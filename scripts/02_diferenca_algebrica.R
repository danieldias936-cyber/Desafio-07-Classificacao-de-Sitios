# =============================================================================
# 02_diferenca_algebrica.R — Parte 2: método da diferença algébrica (ADA)
# -----------------------------------------------------------------------------
# Entrada : dados/Sitio.xls
# Saídas  : resultados/02_diferenca_algebrica/
#   curvas_diferenca_algebrica.csv, indice_local_parcelas_diferenca_algebrica.csv,
#   ajuste_modelos_diferenca_algebrica.csv, coeficientes_diferenca_algebrica.csv,
#   pares_medicoes.csv, feixe_diferenca_algebrica.png, residuos_diferenca_algebrica.png
# Execução: Rscript scripts/02_diferenca_algebrica.R   (a partir da raiz)
# =============================================================================
source("scripts/00_config.R")
saida <- pasta_saida("02_diferenca_algebrica")
d <- le_dados()

# ---- 1. Equações em diferença a partir do Chapman-Richards --------------------
# Modelo de base: Hd = b0 * (1 - exp(-b1*I))^b2
#
# (a) Parâmetro de sítio = b0 (assíntota). Isolando b0 em (I1, H1) e
#     substituindo em I2:
#       H2 = H1 * [(1 - exp(-b1*I2)) / (1 - exp(-b1*I1))]^b2
#     -> curvas ANAMÓRFICAS (assíntotas diferentes, mesmo formato)
#
# (b) Parâmetro de sítio = b2 (forma). Isolando b2 em (I1, H1):
#       H2 = b0 * (H1/b0)^[ ln(1 - exp(-b1*I2)) / ln(1 - exp(-b1*I1)) ]
#     -> curvas POLIMÓRFICAS (assíntota comum, formatos diferentes)
ada_anam <- function(H1, I1, I2, b1, b2) H1 * ((1 - exp(-b1 * I2)) / (1 - exp(-b1 * I1)))^b2
ada_poli <- function(H1, I1, I2, b0, b1) b0 * (H1 / b0)^(log(1 - exp(-b1 * I2)) / log(1 - exp(-b1 * I1)))

# ---- 2. Pares de medições da mesma parcela -----------------------------------
# São usados pares CONSECUTIVOS (1-2, 2-3, ...). Combinar todas as medições
# entre si multiplicaria artificialmente o número de observações e repetiria
# a mesma informação, inflando a precisão aparente do ajuste.
pares <- do.call(rbind, lapply(split(d, d$id), function(x) {
  n <- nrow(x)
  data.frame(id = x$id[1], I1 = x$idade[-n], H1 = x$hd[-n], I2 = x$idade[-1], H2 = x$hd[-1])
}))
cat(sprintf("Pares consecutivos: %d (de %d parcelas)\n", nrow(pares), length(unique(pares$id))))

# ---- 3. Ajuste por regressão não linear --------------------------------------
m_anam <- nls(H2 ~ ada_anam(H1, I1, I2, b1, b2), data = pares,
              start = list(b1 = 0.3, b2 = 1.2), control = nls.control(maxiter = 500))
m_poli <- nls(H2 ~ ada_poli(H1, I1, I2, b0, b1), data = pares,
              start = list(b0 = 35, b1 = 0.3), control = nls.control(maxiter = 500))

tab_ajuste <- rbind(
  metricas(pares$H2, fitted(m_anam), 2, "ADA anamórfica (b0 livre)"),
  metricas(pares$H2, fitted(m_poli), 2, "ADA polimórfica (b2 livre)"))
# Critério adicional: erro de PROJEÇÃO a partir da idade-índice. As curvas de
# sítio são usadas partindo de S (Hd em Ii) para qualquer idade; por isso,
# projeta-se a medição de cada parcela mais próxima de Ii para todas as demais
# idades dessa parcela e calcula-se o RMSE. Pares consecutivos avaliam só
# projeções curtas (~1 ano); este teste avalia projeções de até 3,7 anos.
d_tmp <- d
ref <- medicao_mais_proxima(d_tmp)
proj_teste <- merge(d_tmp, data.frame(id = ref$id, I_ref = ref$idade, H_ref = ref$hd), by = "id")
proj_teste <- proj_teste[proj_teste$idade != proj_teste$I_ref, ]
ca0 <- coef(m_anam); cp0 <- coef(m_poli)
e_anam <- proj_teste$hd - ada_anam(proj_teste$H_ref, proj_teste$I_ref, proj_teste$idade, ca0["b1"], ca0["b2"])
e_poli <- proj_teste$hd - ada_poli(proj_teste$H_ref, proj_teste$I_ref, proj_teste$idade, cp0["b0"], cp0["b1"])
tab_ajuste$rmse_projecao_ii_m <- c(sqrt(mean(e_anam^2)), sqrt(mean(e_poli^2)))
tab_ajuste$vies_projecao_ii_m <- c(mean(e_anam), mean(e_poli))

# Seleção: menor erro de projeção a partir de Ii (uso real das curvas), com o
# Syx dos pares consecutivos como informação complementar.
sel_i <- which.min(tab_ajuste$rmse_projecao_ii_m)
tab_ajuste$selecionado <- seq_len(nrow(tab_ajuste)) == sel_i
print(tab_ajuste, digits = 4)

coefs <- do.call(rbind, lapply(list(m_anam, m_poli), function(m) {
  s <- summary(m)$coefficients
  data.frame(parametro = rownames(s), estimativa = s[, 1], erro_padrao = s[, 2], t = s[, 3], p_valor = s[, 4])
}))
coefs$modelo <- rep(tab_ajuste$modelo, each = 2)

ca <- coef(m_anam); cp <- coef(m_poli)
proj <- if (sel_i == 1) {
  function(H1, I1, I2) ada_anam(H1, I1, I2, ca["b1"], ca["b2"])
} else {
  function(H1, I1, I2) ada_poli(H1, I1, I2, cp["b0"], cp["b1"])
}

# ---- 4. Feixe de curvas: H1 = S em I1 = Ii ------------------------------------
curvas_ambos <- do.call(rbind, lapply(1:2, function(k) {
  g <- expand.grid(idade = IDADES_CURVA, S = CLASSES_S)
  g$hd <- if (k == 1) ada_anam(g$S, IDADE_INDICE, g$idade, ca["b1"], ca["b2"]) else
                      ada_poli(g$S, IDADE_INDICE, g$idade, cp["b0"], cp["b1"])
  g$modelo <- tab_ajuste$modelo[k]
  g
}))
curvas <- subset(curvas_ambos, modelo == tab_ajuste$modelo[sel_i])[, c("idade", "S", "hd")]
curvas$metodo <- "Diferença algébrica"

# ---- 5. Classificação das parcelas --------------------------------------------
d$s_estimado <- as.numeric(proj(d$hd, d$idade, IDADE_INDICE))
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
classif$metodo <- "Diferença algébrica"
print(classif, digits = 4)

# ---- 6. Anamorfismo x polimorfismo: razão entre curvas -----------------------
# Em curvas anamórficas, a razão Hd(maior S)/Hd(menor S) é constante com a idade.
razao <- function(df) {
  a <- df[df$S == max(CLASSES_S), ]; b <- df[df$S == min(CLASSES_S), ]
  range(a$hd / b$hd)
}
cg <- read.csv("resultados/01_curva_guia/curvas_curva_guia.csv")
forma <- data.frame(
  conjunto = c("Curva-guia", tab_ajuste$modelo),
  razao_min = c(razao(cg)[1], razao(subset(curvas_ambos, modelo == tab_ajuste$modelo[1]))[1],
                razao(subset(curvas_ambos, modelo == tab_ajuste$modelo[2]))[1]),
  razao_max = c(razao(cg)[2], razao(subset(curvas_ambos, modelo == tab_ajuste$modelo[1]))[2],
                razao(subset(curvas_ambos, modelo == tab_ajuste$modelo[2]))[2]))
forma$tipo <- ifelse(forma$razao_max - forma$razao_min < 1e-6, "anamórfico", "polimórfico")
print(forma, digits = 4)

# ---- 7. Saídas ----------------------------------------------------------------
write.csv(curvas,     file.path(saida, "curvas_diferenca_algebrica.csv"), row.names = FALSE)
write.csv(classif,    file.path(saida, "indice_local_parcelas_diferenca_algebrica.csv"), row.names = FALSE)
write.csv(tab_ajuste, file.path(saida, "ajuste_modelos_diferenca_algebrica.csv"), row.names = FALSE)
write.csv(coefs,      file.path(saida, "coeficientes_diferenca_algebrica.csv"), row.names = FALSE)
write.csv(pares,      file.path(saida, "pares_medicoes.csv"), row.names = FALSE)
write.csv(forma,      file.path(saida, "anamorfismo_polimorfismo.csv"), row.names = FALSE)

res <- rbind(data.frame(modelo = tab_ajuste$modelo[1], estimado = fitted(m_anam), residuo = resid(m_anam)),
             data.frame(modelo = tab_ajuste$modelo[2], estimado = fitted(m_poli), residuo = resid(m_poli)))
ggsave(file.path(saida, "residuos_diferenca_algebrica.png"),
       grafico_residuos(res, "Diferença algébrica: resíduos (Hd2 observada - estimada)"),
       width = 9, height = 3.8, dpi = 200, bg = "white")
ggsave(file.path(saida, "feixe_diferenca_algebrica.png"),
       grafico_feixe(curvas, d, sprintf("Diferença algébrica — %s", tab_ajuste$modelo[sel_i])),
       width = 8, height = 5.5, dpi = 200, bg = "white")
g2 <- ggplot(curvas_ambos, aes(idade, hd, group = factor(S))) + geom_line(colour = "grey30") +
  facet_wrap(~ modelo) + geom_point(data = d, aes(idade, hd), inherit.aes = FALSE, colour = "#2E7D32", size = 1.5) +
  labs(x = "Idade (anos)", y = "Altura dominante (m)", title = "Formulações anamórfica e polimórfica") + tema
ggsave(file.path(saida, "comparacao_formulacoes.png"), g2, width = 9, height = 4.2, dpi = 200, bg = "white")
cat("Concluído: resultados/02_diferenca_algebrica/\n")
