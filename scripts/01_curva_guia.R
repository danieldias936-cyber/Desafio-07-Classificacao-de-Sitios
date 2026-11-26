# =============================================================================
# 01_curva_guia.R — Parte 1: método da curva-guia
# -----------------------------------------------------------------------------
# Entrada : dados/Sitio.xls
# Saídas  : resultados/01_curva_guia/
#   curvas_curva_guia.csv, indice_local_parcelas_curva_guia.csv,
#   ajuste_modelos_curva_guia.csv, coeficientes_curva_guia.csv,
#   feixe_curva_guia.png, residuos_curva_guia.png
# Execução: Rscript scripts/01_curva_guia.R   (a partir da raiz)
# =============================================================================
source("scripts/00_config.R")
saida <- pasta_saida("01_curva_guia")
d <- le_dados()

cat(sprintf("Base: %d medições, %d parcelas, idade %.2f-%.2f anos, Hd %.1f-%.1f m\n",
            nrow(d), length(unique(d$id)), min(d$idade), max(d$idade), min(d$hd), max(d$hd)))

# ---- 1. Ajuste de três modelos com todos os dados agrupados -----------------
modelos <- list(
  "Chapman-Richards" = list(
    f = hd ~ b0 * (1 - exp(-b1 * idade))^b2,
    ini = list(b0 = 30, b1 = 0.3, b2 = 1.2),
    pred = function(cf, I) cf["b0"] * (1 - exp(-cf["b1"] * I))^cf["b2"]),
  "Schumacher" = list(
    f = hd ~ b0 * exp(b1 / idade),
    ini = list(b0 = 30, b1 = -2),
    pred = function(cf, I) cf["b0"] * exp(cf["b1"] / I)),
  "Logístico" = list(
    f = hd ~ b0 / (1 + b1 * exp(-b2 * idade)),
    ini = list(b0 = 25, b1 = 3, b2 = 0.6),
    pred = function(cf, I) cf["b0"] / (1 + cf["b1"] * exp(-cf["b2"] * I)))
)

ajustes <- lapply(modelos, function(m) nls(m$f, data = d, start = m$ini,
                                          control = nls.control(maxiter = 500)))
tab_ajuste <- do.call(rbind, lapply(names(ajustes), function(nm) {
  cf <- coef(ajustes[[nm]])
  metricas(d$hd, modelos[[nm]]$pred(cf, d$idade), length(cf), nm)
}))
print(tab_ajuste, digits = 4)

# ---- 2. Seleção: Chapman-Richards como modelo de base ------------------------
# O enunciado exige o MESMO modelo de base nas partes 1, 2, 3 e 5. O
# Chapman-Richards é mantido quando seu desempenho é equivalente ao do melhor
# modelo (diferença de Syx% menor que 0,5 ponto percentual); caso contrário,
# o melhor modelo seria usado na curva-guia e a diferença, documentada.
melhor <- tab_ajuste$modelo[which.min(tab_ajuste$syx_pct)]
dif_cr <- tab_ajuste$syx_pct[tab_ajuste$modelo == "Chapman-Richards"] - min(tab_ajuste$syx_pct)
sel <- if (dif_cr < 0.5) "Chapman-Richards" else melhor
tab_ajuste$selecionado <- tab_ajuste$modelo == sel
cat(sprintf("Melhor Syx%%: %s | diferença do Chapman-Richards: %.2f p.p. | selecionado: %s\n",
            melhor, dif_cr, sel))

cf  <- coef(ajustes[[sel]])
f_g <- function(I) as.numeric(modelos[[sel]]$pred(cf, I))   # curva-guia

coefs <- as.data.frame(summary(ajustes[[sel]])$coefficients)
coefs <- data.frame(modelo = sel, parametro = rownames(coefs), estimativa = coefs[, 1],
                    erro_padrao = coefs[, 2], t = coefs[, 3], p_valor = coefs[, 4])

# ---- 3. Relação entre S e Hd (método da proporcionalidade) -------------------
#   Hd(I) = S * f(I) / f(Ii)      e      S = Hd * f(Ii) / f(I)
s_de_hd <- function(hd, I) hd * f_g(IDADE_INDICE) / f_g(I)
hd_de_s <- function(S, I)  S * f_g(I) / f_g(IDADE_INDICE)

curvas <- expand.grid(idade = IDADES_CURVA, S = CLASSES_S)
curvas$hd <- hd_de_s(curvas$S, curvas$idade)
curvas$metodo <- "Curva-guia"

# ---- 4. Classificação das 9 parcelas ------------------------------------------
d$s_estimado <- s_de_hd(d$hd, d$idade)
prox <- medicao_mais_proxima(d)
estab <- do.call(rbind, lapply(split(d, d$id), function(x) {
  cl <- classe_s(x$s_estimado)
  data.frame(id = x$id[1], s_min = min(x$s_estimado), s_max = max(x$s_estimado),
             trocas_de_classe = sum(diff(cl) != 0))
}))
classif <- merge(data.frame(id = prox$id, talhao = prox$talhao, parcela = prox$parcela,
                            idade_usada = prox$idade, hd_usada = prox$hd,
                            S = prox$s_estimado, classe_S = classe_s(prox$s_estimado)),
                 estab, by = "id")
classif$metodo <- "Curva-guia"
print(classif, digits = 4)

# ---- 5. Saídas ----------------------------------------------------------------
write.csv(curvas,     file.path(saida, "curvas_curva_guia.csv"), row.names = FALSE)
write.csv(classif,    file.path(saida, "indice_local_parcelas_curva_guia.csv"), row.names = FALSE)
write.csv(tab_ajuste, file.path(saida, "ajuste_modelos_curva_guia.csv"), row.names = FALSE)
write.csv(coefs,      file.path(saida, "coeficientes_curva_guia.csv"), row.names = FALSE)

res <- do.call(rbind, lapply(names(ajustes), function(nm) {
  est <- as.numeric(modelos[[nm]]$pred(coef(ajustes[[nm]]), d$idade))
  data.frame(modelo = nm, estimado = est, residuo = d$hd - est)
}))
ggsave(file.path(saida, "residuos_curva_guia.png"),
       grafico_residuos(res, "Curva-guia: resíduos dos modelos testados"),
       width = 10, height = 3.8, dpi = 200, bg = "white")
ggsave(file.path(saida, "feixe_curva_guia.png"),
       grafico_feixe(curvas, d, sprintf("Curva-guia (%s)", sel)),
       width = 8, height = 5.5, dpi = 200, bg = "white")

# Gráfico exploratório inicial (Hd x idade por parcela)
g0 <- ggplot(d, aes(idade, hd, colour = id)) + geom_line() + geom_point(size = 2) +
  labs(x = "Idade (anos)", y = "Altura dominante (m)", colour = "Parcela",
       title = "Trajetórias observadas de altura dominante") + tema
ggsave(file.path(saida, "exploratorio_hd_idade.png"), g0, width = 8, height = 5, dpi = 200, bg = "white")

cat("Concluído: resultados/01_curva_guia/\n")
