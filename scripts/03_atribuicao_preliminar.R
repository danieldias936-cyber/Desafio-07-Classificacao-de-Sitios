# =============================================================================
# 03_atribuicao_preliminar.R — Parte 3: atribuição preliminar de índices de local
# -----------------------------------------------------------------------------
# Entradas: dados/Sitio.xls
#           resultados/01_curva_guia/indice_local_parcelas_curva_guia.csv
# Saídas  : resultados/03_atribuicao_preliminar/
#   curvas_atribuicao_preliminar.csv, indice_local_parcelas_atribuicao_preliminar.csv,
#   ajuste_modelos_atribuicao_preliminar.csv, coeficientes_atribuicao_preliminar.csv,
#   feixe_atribuicao_preliminar.png, residuos_atribuicao_preliminar.png
# Execução: Rscript scripts/03_atribuicao_preliminar.R   (após o script 01)
# =============================================================================
source("scripts/00_config.R")
saida <- pasta_saida("03_atribuicao_preliminar")
d <- le_dados()

# ---- 1. Índice de local preliminar (curva-guia, script 01) --------------------
arq01 <- "resultados/01_curva_guia/indice_local_parcelas_curva_guia.csv"
if (!file.exists(arq01)) stop("Rode antes o script 01_curva_guia.R")
prelim <- read.csv(arq01)[, c("id", "S")]
names(prelim)[2] <- "S_prelim"
d <- merge(d, prelim, by = "id")
d <- d[order(d$id, d$idade), ]

# ---- 2. Modelo de Hd em função da idade e de S --------------------------------
# Forma do Chapman-Richards reescrita para passar por S na idade-índice:
#   Hd = S * [ (1 - exp(-k*I)) / (1 - exp(-k*Ii)) ]^b2
# (a) anamórfico:  k = c0                 (S só desloca a curva)
# (b) polimórfico: k = c0 + c1*S          (S altera também a forma)
# Por construção, Hd(Ii) = S em ambos, o que mantém o significado de S.
m_ap <- function(I, S, c0, c1, b2) {
  k <- c0 + c1 * S
  S * ((1 - exp(-k * I)) / (1 - exp(-k * IDADE_INDICE)))^b2
}
fit_a <- nls(hd ~ m_ap(idade, S_prelim, c0, 0, b2), data = d,
             start = list(c0 = 0.3, b2 = 1.2), control = nls.control(maxiter = 500))
fit_p <- nls(hd ~ m_ap(idade, S_prelim, c0, c1, b2), data = d,
             start = list(c0 = 0.3, c1 = 0, b2 = 1.2), control = nls.control(maxiter = 500))

tab_ajuste <- rbind(metricas(d$hd, fitted(fit_a), 2, "Anamórfico (k constante)"),
                    metricas(d$hd, fitted(fit_p), 3, "Polimórfico (k = c0 + c1*S)"))
teste_f <- anova(fit_a, fit_p)
tab_ajuste$p_valor_teste_F <- c(NA, teste_f$`Pr(>F)`[2])
# Seleção: o modelo polimórfico só é adotado se o termo c1 for significativo
# (teste F de modelos aninhados, 5%); caso contrário, fica o mais simples.
sel_i <- if (!is.na(teste_f$`Pr(>F)`[2]) && teste_f$`Pr(>F)`[2] < 0.05) 2 else 1
tab_ajuste$selecionado <- seq_len(2) == sel_i
print(tab_ajuste, digits = 4)

coefs <- do.call(rbind, lapply(list(fit_a, fit_p), function(m) {
  s <- summary(m)$coefficients
  data.frame(parametro = rownames(s), estimativa = s[, 1], erro_padrao = s[, 2], t = s[, 3], p_valor = s[, 4])
}))
coefs$modelo <- c(rep(tab_ajuste$modelo[1], 2), rep(tab_ajuste$modelo[2], 3))

cf <- coef(if (sel_i == 1) fit_a else fit_p)
c1 <- if (sel_i == 1) 0 else cf["c1"]
hd_fun <- function(I, S) as.numeric(m_ap(I, S, cf["c0"], c1, cf["b2"]))

# ---- 3. Feixe de curvas --------------------------------------------------------
curvas <- expand.grid(idade = IDADES_CURVA, S = CLASSES_S)
curvas$hd <- hd_fun(curvas$idade, curvas$S)
curvas$metodo <- "Atribuição preliminar"

# ---- 4. Reclassificação com o modelo final ------------------------------------
# Para cada medição, resolve-se Hd_obs = f(I, S) em S (uniroot).
s_inv <- function(hd, I) uniroot(function(S) hd_fun(I, S) - hd, c(5, 60))$root
d$s_estimado <- mapply(s_inv, d$hd, d$idade)
prox <- medicao_mais_proxima(d)
estab <- do.call(rbind, lapply(split(d, d$id), function(x) {
  cl <- classe_s(x$s_estimado)
  data.frame(id = x$id[1], s_min = min(x$s_estimado), s_max = max(x$s_estimado),
             trocas_de_classe = sum(diff(cl) != 0))
}))
classif <- merge(data.frame(id = prox$id, talhao = prox$talhao, parcela = prox$parcela,
                            idade_usada = prox$idade, hd_usada = prox$hd,
                            S_preliminar = prox$S_prelim, classe_S_preliminar = classe_s(prox$S_prelim),
                            S = prox$s_estimado, classe_S = classe_s(prox$s_estimado)),
                 estab, by = "id")
classif$diferenca_S <- classif$S - classif$S_preliminar
classif$mudou_classe <- classif$classe_S != classif$classe_S_preliminar
classif$metodo <- "Atribuição preliminar"
print(classif[, c("id", "S_preliminar", "S", "diferenca_S", "classe_S_preliminar", "classe_S", "mudou_classe")], digits = 4)
cat(sprintf("Parcelas que mudaram de classe: %d de %d\n", sum(classif$mudou_classe), nrow(classif)))

# ---- 5. Saídas ----------------------------------------------------------------
write.csv(curvas,     file.path(saida, "curvas_atribuicao_preliminar.csv"), row.names = FALSE)
write.csv(classif,    file.path(saida, "indice_local_parcelas_atribuicao_preliminar.csv"), row.names = FALSE)
write.csv(tab_ajuste, file.path(saida, "ajuste_modelos_atribuicao_preliminar.csv"), row.names = FALSE)
write.csv(coefs,      file.path(saida, "coeficientes_atribuicao_preliminar.csv"), row.names = FALSE)

res <- rbind(data.frame(modelo = tab_ajuste$modelo[1], estimado = fitted(fit_a), residuo = resid(fit_a)),
             data.frame(modelo = tab_ajuste$modelo[2], estimado = fitted(fit_p), residuo = resid(fit_p)))
ggsave(file.path(saida, "residuos_atribuicao_preliminar.png"),
       grafico_residuos(res, "Atribuição preliminar: resíduos"),
       width = 9, height = 3.8, dpi = 200, bg = "white")
ggsave(file.path(saida, "feixe_atribuicao_preliminar.png"),
       grafico_feixe(curvas, d, sprintf("Atribuição preliminar — %s", tab_ajuste$modelo[sel_i])),
       width = 8, height = 5.5, dpi = 200, bg = "white")
cat("Concluído: resultados/03_atribuicao_preliminar/\n")
