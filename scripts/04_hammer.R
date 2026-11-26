# =============================================================================
# 04_hammer.R — Parte 4: método de Hammer
# -----------------------------------------------------------------------------
# Entrada : dados/Sitio.xls
# Saídas  : resultados/04_hammer/
#   curvas_hammer.csv, indice_local_parcelas_hammer.csv, pontos_representativos_hammer.csv,
#   ajuste_modelos_hammer.csv, coeficientes_hammer.csv,
#   feixe_hammer.png, pontos_hammer.png, residuos_hammer.png
# Execução: Rscript scripts/04_hammer.R   (a partir da raiz)
# =============================================================================
source("scripts/00_config.R")
saida <- pasta_saida("04_hammer")
d <- le_dados()

# ---- 1. Classes de idade (não por parcela) -----------------------------------
# As medições ocorreram em cinco ocasiões; as classes envolvem cada ocasião.
limites <- c(2.5, 3.0, 4.0, 5.0, 6.0, 6.5)
rotulos <- c("2,5-3,0", "3,5-4,0", "4,5-5,0", "5,5-6,0", "6,0-6,5")
d$classe_idade <- cut(d$idade, breaks = limites, labels = rotulos, right = FALSE)
print(table(d$classe_idade))

# ---- 2. Classes de produtividade dentro de cada classe de idade --------------
# Em cada classe de idade, as medições são ordenadas por Hd e divididas em
# terços (inferior, médio, superior). A Hd representativa de cada classe é a
# média da classe; a idade representativa, a média das idades.
d <- do.call(rbind, lapply(split(d, d$classe_idade), function(x) {
  r <- rank(x$hd, ties.method = "first")
  x$classe_prod <- factor(c("Inferior", "Médio", "Superior")[ceiling(3 * r / nrow(x))],
                          levels = c("Inferior", "Médio", "Superior"))
  x
}))
pontos <- aggregate(cbind(idade, hd) ~ classe_idade + classe_prod, data = d, FUN = mean)
pontos$n <- aggregate(hd ~ classe_idade + classe_prod, data = d, FUN = length)$hd
pontos <- pontos[order(pontos$classe_prod, pontos$idade), ]
print(pontos, digits = 4)

# ---- 3. Ajuste das curvas por classe de produtividade ------------------------
# (a) Chapman-Richards com assíntota própria por classe e forma comum
#     (estratificação do método, curvas anamórficas): 5 parâmetros, 15 pontos
# (b) Chapman-Richards ajustado separadamente em cada classe (3 parâmetros
#     com 5 pontos cada), testado para comparação
pontos$cp <- as.integer(pontos$classe_prod)
fit_comum <- nls(hd ~ b0[cp] * (1 - exp(-b1 * idade))^b2, data = pontos,
                 start = list(b0 = c(25, 28, 32), b1 = 0.3, b2 = 1.2),
                 control = nls.control(maxiter = 500))
fit_sep <- lapply(levels(pontos$classe_prod), function(k) {
  tryCatch(nls(hd ~ b0 * (1 - exp(-b1 * idade))^b2, data = subset(pontos, classe_prod == k),
               start = list(b0 = 28, b1 = 0.3, b2 = 1.2), control = nls.control(maxiter = 500)),
           error = function(e) NULL)
})
names(fit_sep) <- levels(pontos$classe_prod)
sep_ok <- all(!sapply(fit_sep, is.null))
est_sep <- if (sep_ok) unlist(lapply(levels(pontos$classe_prod), function(k)
  fitted(fit_sep[[k]]))) else rep(NA, nrow(pontos))

tab_ajuste <- rbind(metricas(pontos$hd, fitted(fit_comum), 5, "Assíntota por classe, forma comum"),
                    if (sep_ok) metricas(pontos$hd, est_sep, 9, "Curva separada por classe"))
# Seleção: as curvas separadas usam 9 parâmetros para 15 pontos (3 para 5
# pontos em cada classe). Só são preferidas se reduzirem o AIC em mais de 2
# unidades; diferenças menores indicam modelos equivalentes, e fica o mais
# parcimonioso (forma comum).
sel_i <- if (sep_ok && (tab_ajuste$aic[1] - tab_ajuste$aic[2]) > 2) 2 else 1
tab_ajuste$selecionado <- seq_len(nrow(tab_ajuste)) == sel_i
print(tab_ajuste, digits = 4)
if (!sep_ok) cat("Aviso: o ajuste separado não convergiu em todas as classes.\n")

coefs <- as.data.frame(summary(fit_comum)$coefficients)
coefs <- data.frame(modelo = tab_ajuste$modelo[1], parametro = c("b0_inferior", "b0_medio", "b0_superior", "b1", "b2"),
                    estimativa = coefs[, 1], erro_padrao = coefs[, 2], t = coefs[, 3], p_valor = coefs[, 4])

# O feixe é gerado com a forma comum (b1, b2) do modelo selecionado
cf <- coef(fit_comum)
forma <- function(I) (1 - exp(-cf["b1"] * I))^cf["b2"]
s_classes <- data.frame(classe_prod = levels(pontos$classe_prod),
                        b0 = cf[c("b01", "b02", "b03")],
                        S_na_idade_indice = cf[c("b01", "b02", "b03")] * forma(IDADE_INDICE))
print(s_classes, digits = 4)

# ---- 4. Feixe de curvas para as classes de S ----------------------------------
curvas <- expand.grid(idade = IDADES_CURVA, S = CLASSES_S)
curvas$hd <- as.numeric(curvas$S * forma(curvas$idade) / forma(IDADE_INDICE))
curvas$metodo <- "Hammer"

# ---- 5. Classificação das parcelas --------------------------------------------
d$s_estimado <- as.numeric(d$hd * forma(IDADE_INDICE) / forma(d$idade))
prox <- medicao_mais_proxima(d)
estab <- do.call(rbind, lapply(split(d, d$id), function(x) {
  cl <- classe_s(x$s_estimado)
  tp <- table(x$classe_prod)
  data.frame(id = x$id[1], s_min = min(x$s_estimado), s_max = max(x$s_estimado),
             trocas_de_classe = sum(diff(cl) != 0),
             classe_prod_predominante = names(tp)[which.max(tp)])
}))
classif <- merge(data.frame(id = prox$id, talhao = prox$talhao, parcela = prox$parcela,
                            idade_usada = prox$idade, hd_usada = prox$hd,
                            S = prox$s_estimado, classe_S = classe_s(prox$s_estimado)),
                 estab, by = "id")
classif$metodo <- "Hammer"
print(classif, digits = 4)

# ---- 6. Saídas ----------------------------------------------------------------
write.csv(curvas,     file.path(saida, "curvas_hammer.csv"), row.names = FALSE)
write.csv(classif,    file.path(saida, "indice_local_parcelas_hammer.csv"), row.names = FALSE)
write.csv(pontos[, c("classe_idade", "classe_prod", "n", "idade", "hd")],
          file.path(saida, "pontos_representativos_hammer.csv"), row.names = FALSE)
write.csv(tab_ajuste, file.path(saida, "ajuste_modelos_hammer.csv"), row.names = FALSE)
write.csv(coefs,      file.path(saida, "coeficientes_hammer.csv"), row.names = FALSE)

res <- data.frame(modelo = tab_ajuste$modelo[1], estimado = fitted(fit_comum), residuo = resid(fit_comum))
if (sep_ok) res <- rbind(res, data.frame(modelo = tab_ajuste$modelo[2], estimado = est_sep, residuo = pontos$hd - est_sep))
ggsave(file.path(saida, "residuos_hammer.png"), grafico_residuos(res, "Hammer: resíduos nos pontos representativos"),
       width = 9, height = 3.8, dpi = 200, bg = "white")
ggsave(file.path(saida, "feixe_hammer.png"),
       grafico_feixe(curvas, d, "Método de Hammer (forma comum às classes de produtividade)"),
       width = 8, height = 5.5, dpi = 200, bg = "white")
grade <- expand.grid(idade = IDADES_CURVA, cp = 1:3)
grade$hd <- cf[c("b01", "b02", "b03")][grade$cp] * forma(grade$idade)
grade$classe_prod <- levels(pontos$classe_prod)[grade$cp]
gp <- ggplot() +
  geom_point(data = d, aes(idade, hd, colour = classe_prod), alpha = 0.35) +
  geom_point(data = pontos, aes(idade, hd, colour = classe_prod), size = 3.5, shape = 17) +
  geom_line(data = grade, aes(idade, hd, colour = classe_prod), linewidth = 1) +
  labs(x = "Idade (anos)", y = "Altura dominante (m)", colour = "Classe de\nprodutividade",
       title = "Hammer: pontos representativos (triângulos) e curvas por classe") + tema
ggsave(file.path(saida, "pontos_hammer.png"), gp, width = 8, height = 5.5, dpi = 200, bg = "white")
cat("Concluído: resultados/04_hammer/\n")
