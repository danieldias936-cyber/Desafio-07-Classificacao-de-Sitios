# =============================================================================
# 05_predicao_parametros.R — Parte 5: método da predição de parâmetros
# -----------------------------------------------------------------------------
# Entrada : dados/Sitio.xls
# Saídas  : resultados/05_predicao_parametros/
#   curvas_predicao_parametros.csv, indice_local_parcelas_predicao_parametros.csv,
#   parametros_por_parcela.csv, relacoes_parametros_S.csv, ajuste_modelos_predicao_parametros.csv,
#   feixe_predicao_parametros.png, parametros_vs_S.png, residuos_predicao_parametros.png
# Execução: Rscript scripts/05_predicao_parametros.R   (a partir da raiz)
# =============================================================================
source("scripts/00_config.R")
saida <- pasta_saida("05_predicao_parametros")
d <- le_dados()

# ---- 1. Ajuste do Chapman-Richards por parcela --------------------------------
# Cada parcela tem só 4 ou 5 medições. O modelo completo (3 parâmetros por
# parcela) é testado, mas gera assíntotas instáveis (b0 de 21 a 41 m) porque
# b0, b1 e b2 se compensam com tão poucos pontos. Por isso, um parâmetro é
# fixado no valor do ajuste conjunto e os outros dois são estimados por
# parcela. Duas alternativas são comparadas:
#   A: b2 fixo; b0 e b1 por parcela
#   B: b1 fixo; b0 e b2 por parcela
global <- nls(hd ~ b0 * (1 - exp(-b1 * idade))^b2, data = d,
              start = list(b0 = 30, b1 = 0.3, b2 = 1.2), control = nls.control(maxiter = 500))
G <- coef(global)

completo <- sapply(split(d, d$id), function(x) {
  f <- tryCatch(nls(hd ~ b0 * (1 - exp(-b1 * idade))^b2, data = x,
                    start = list(b0 = 30, b1 = 0.3, b2 = 1.2), control = nls.control(maxiter = 500)),
                error = function(e) NULL)
  if (is.null(f)) NA else unname(coef(f)["b0"])
})

ajusta_parcelas <- function(livres) {
  do.call(rbind, lapply(split(d, d$id), function(x) {
    if (livres == "b0_b1") {
      f <- nls(hd ~ b0 * (1 - exp(-b1 * idade))^G[["b2"]], data = x,
               start = list(b0 = 30, b1 = 0.3), control = nls.control(maxiter = 500))
      cf <- c(coef(f), b2 = G[["b2"]])
    } else {
      f <- nls(hd ~ b0 * (1 - exp(-G[["b1"]] * idade))^b2, data = x,
               start = list(b0 = 30, b2 = 1), control = nls.control(maxiter = 500))
      cf <- c(b0 = coef(f)[["b0"]], b1 = G[["b1"]], b2 = coef(f)[["b2"]])
    }
    data.frame(id = x$id[1], talhao = x$talhao[1], parcela = x$parcela[1], n = nrow(x),
               b0 = cf[["b0"]], b1 = cf[["b1"]], b2 = cf[["b2"]],
               syx_m = sqrt(sum(resid(f)^2) / (nrow(x) - 2)),
               S = chapman(IDADE_INDICE, cf[["b0"]], cf[["b1"]], cf[["b2"]]))
  }))
}

# ---- 2. Relações parâmetro ~ S e equação única Hd = f(I, S) -------------------
# Cada parâmetro livre é regredido linearmente em S. Se a relação não for
# significativa (5%), o parâmetro é tratado como constante (média).
monta <- function(livres) {
  pp <- ajusta_parcelas(livres)
  vars <- if (livres == "b0_b1") c("b0", "b1") else c("b0", "b2")
  rel <- do.call(rbind, lapply(vars, function(v) {
    m <- lm(as.formula(paste(v, "~ S")), data = pp)
    s <- coef(summary(m))
    data.frame(alternativa = livres, parametro = v, intercepto = s[1, 1], inclinacao = s[2, 1],
               p_valor_inclinacao = s[2, 4], r2 = summary(m)$r.squared, media = mean(pp[[v]]))
  }))
  rel$usa_S <- rel$p_valor_inclinacao < 0.05
  par_S <- function(v, S) {
    r <- rel[rel$parametro == v, ]
    if (nrow(r) == 0) return(G[[v]])
    if (r$usa_S) r$intercepto + r$inclinacao * S else r$media
  }
  f <- function(I, S) as.numeric(par_S("b0", S) * (1 - exp(-par_S("b1", S) * I))^par_S("b2", S))
  est <- f(d$idade, pp$S[match(d$id, pp$id)])
  met <- metricas(d$hd, est, 2 + sum(rel$usa_S) * 2 + sum(!rel$usa_S), paste("Alternativa", ifelse(livres == "b0_b1", "A (b0, b1 livres)", "B (b0, b2 livres)")))
  met$desvio_max_em_Ii_m <- max(abs(f(IDADE_INDICE, CLASSES_S) - CLASSES_S))
  met$relacoes_significativas <- sum(rel$usa_S)
  list(pp = pp, rel = rel, f = f, met = met, est = est)
}
alt <- list(A = monta("b0_b1"), B = monta("b0_b2"))
tab_ajuste <- rbind(alt$A$met, alt$B$met)

# Seleção: a equação final precisa reproduzir S na idade-índice (por
# definição, Hd(Ii) = S); escolhe-se a alternativa com menor desvio em Ii,
# desde que tenha ao menos uma relação significativa com S.
ok <- tab_ajuste$relacoes_significativas > 0
sel <- if (any(ok)) which(ok)[which.min(tab_ajuste$desvio_max_em_Ii_m[ok])] else which.min(tab_ajuste$desvio_max_em_Ii_m)
tab_ajuste$selecionado <- seq_len(2) == sel
print(tab_ajuste, digits = 4)
E <- alt[[sel]]
rel <- rbind(alt$A$rel, alt$B$rel)
print(rel, digits = 4)

par_parc <- E$pp
par_parc$b0_modelo_completo <- completo[par_parc$id]
print(par_parc, digits = 4)

# ---- 3. Feixe de curvas --------------------------------------------------------
curvas <- expand.grid(idade = IDADES_CURVA, S = CLASSES_S)
curvas$hd <- E$f(curvas$idade, curvas$S)
curvas$metodo <- "Predição de parâmetros"

# ---- 4. Classificação das parcelas --------------------------------------------
s_inv <- function(hd, I) uniroot(function(S) E$f(I, S) - hd, c(5, 60))$root
d$s_estimado <- mapply(s_inv, d$hd, d$idade)
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
classif$S_curva_propria <- par_parc$S[match(classif$id, par_parc$id)]
classif$metodo <- "Predição de parâmetros"
print(classif, digits = 4)

# ---- 5. Saídas ----------------------------------------------------------------
write.csv(curvas,     file.path(saida, "curvas_predicao_parametros.csv"), row.names = FALSE)
write.csv(classif,    file.path(saida, "indice_local_parcelas_predicao_parametros.csv"), row.names = FALSE)
write.csv(par_parc,   file.path(saida, "parametros_por_parcela.csv"), row.names = FALSE)
write.csv(rel,        file.path(saida, "relacoes_parametros_S.csv"), row.names = FALSE)
write.csv(tab_ajuste, file.path(saida, "ajuste_modelos_predicao_parametros.csv"), row.names = FALSE)

res <- rbind(data.frame(modelo = tab_ajuste$modelo[1], estimado = alt$A$est, residuo = d$hd - alt$A$est),
             data.frame(modelo = tab_ajuste$modelo[2], estimado = alt$B$est, residuo = d$hd - alt$B$est))
ggsave(file.path(saida, "residuos_predicao_parametros.png"),
       grafico_residuos(res, "Predição de parâmetros: resíduos da equação única"),
       width = 9, height = 3.8, dpi = 200, bg = "white")
lp <- do.call(rbind, lapply(c("b0", "b1", "b2"), function(v)
  if (sd(par_parc[[v]]) > 0) data.frame(parametro = v, S = par_parc$S, valor = par_parc[[v]])))
gpar <- ggplot(lp, aes(S, valor)) + geom_point(size = 2.5, colour = "#2E7D32") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey30") +
  facet_wrap(~ parametro, scales = "free_y") +
  labs(x = "Índice de local da parcela, S (m)", y = "Estimativa por parcela",
       title = "Parâmetros estimados por parcela em função de S") + tema
ggsave(file.path(saida, "parametros_vs_S.png"), gpar, width = 9, height = 4, dpi = 200, bg = "white")
ggsave(file.path(saida, "feixe_predicao_parametros.png"),
       grafico_feixe(curvas, d, sprintf("Predição de parâmetros — %s", tab_ajuste$modelo[sel])),
       width = 8, height = 5.5, dpi = 200, bg = "white")
cat("Concluído: resultados/05_predicao_parametros/\n")
