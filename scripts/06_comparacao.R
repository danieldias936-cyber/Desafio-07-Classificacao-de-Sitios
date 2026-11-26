# =============================================================================
# 06_comparacao.R — consolidação dos cinco métodos (para o relatório)
# -----------------------------------------------------------------------------
# Entradas: saídas dos scripts 01 a 05
# Saídas  : resultados/06_comparacao/
#   comparacao_curvas.png, comparacao_classificacao.csv, resumo_metodos.csv
# Execução: Rscript scripts/06_comparacao.R   (após os scripts 01 a 05)
# =============================================================================
source("scripts/00_config.R")
saida <- pasta_saida("06_comparacao")
d <- le_dados()

met <- c("01_curva_guia" = "curva_guia", "02_diferenca_algebrica" = "diferenca_algebrica",
         "03_atribuicao_preliminar" = "atribuicao_preliminar", "04_hammer" = "hammer",
         "05_predicao_parametros" = "predicao_parametros")
for (p in names(met)) if (!file.exists(file.path("resultados", p, paste0("curvas_", met[p], ".csv"))))
  stop("Rode antes os scripts 01 a 05.")

curvas <- do.call(rbind, lapply(names(met), function(p)
  read.csv(file.path("resultados", p, paste0("curvas_", met[p], ".csv")))[, c("idade", "S", "hd", "metodo")]))
ordem <- c("Curva-guia", "Diferença algébrica", "Atribuição preliminar", "Hammer", "Predição de parâmetros")
curvas$metodo <- factor(curvas$metodo, levels = ordem)

g <- ggplot(curvas, aes(idade, hd, colour = metodo, group = interaction(metodo, S))) +
  geom_line(linewidth = 0.8) +
  geom_point(data = d, aes(idade, hd), inherit.aes = FALSE, colour = "grey55", size = 1.6) +
  geom_vline(xintercept = IDADE_INDICE, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 6.55, y = CLASSES_S, label = paste0("S = ", CLASSES_S), hjust = 0, size = 3.2) +
  scale_colour_manual(values = c("#1B5E20", "#1565C0", "#EF6C00", "#6A1B9A", "#C62828")) +
  scale_x_continuous(limits = c(2.4, 7.1), breaks = 2:7) +
  labs(x = "Idade (anos)", y = "Altura dominante (m)", colour = "Método",
       title = "Feixes de curvas de índice de local pelos cinco métodos",
       subtitle = sprintf("Idade-índice = %g anos; pontos = medições das 9 parcelas", IDADE_INDICE)) +
  tema + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(saida, "comparacao_curvas.png"), g, width = 9, height = 6.5, dpi = 200, bg = "white")

# Classificação das parcelas pelos cinco métodos
cl <- do.call(rbind, lapply(names(met), function(p) {
  x <- read.csv(file.path("resultados", p, paste0("indice_local_parcelas_", met[p], ".csv")))
  data.frame(id = x$id, metodo = x$metodo, S = x$S, classe_S = x$classe_S, trocas = x$trocas_de_classe)
}))
largo_S  <- reshape(cl[, c("id", "metodo", "S")], idvar = "id", timevar = "metodo", direction = "wide")
largo_cl <- reshape(cl[, c("id", "metodo", "classe_S")], idvar = "id", timevar = "metodo", direction = "wide")
names(largo_S)  <- sub("^S\\.", "S | ", names(largo_S))
names(largo_cl) <- sub("^classe_S\\.", "classe | ", names(largo_cl))
comp <- merge(largo_S, largo_cl, by = "id")
comp$classes_distintas <- apply(largo_cl[match(comp$id, largo_cl$id), -1], 1, function(v) length(unique(v)))
comp$amplitude_S <- apply(largo_S[match(comp$id, largo_S$id), -1], 1, function(v) diff(range(v)))
write.csv(comp, file.path(saida, "comparacao_classificacao.csv"), row.names = FALSE)
print(comp[, c("id", grep("^classe", names(comp), value = TRUE), "amplitude_S")])

# Resumo: estabilidade (trocas de classe somadas) e métricas do modelo selecionado
ajuste <- function(p, arq) {
  a <- read.csv(file.path("resultados", p, arq))
  if ("selecionado" %in% names(a)) a <- a[a$selecionado, ]
  a[1, c("modelo", "n", "r2", "syx_m", "syx_pct")]
}
resumo <- do.call(rbind, Map(function(p, nm) {
  a <- ajuste(p, paste0("ajuste_modelos_", met[p], ".csv"))
  data.frame(metodo = nm, modelo_selecionado = a$modelo, n_ajuste = a$n,
             r2 = a$r2, syx_m = a$syx_m, syx_pct = a$syx_pct,
             trocas_de_classe_total = sum(cl$trocas[cl$metodo == nm]))
}, names(met), ordem))
rownames(resumo) <- NULL
write.csv(resumo, file.path(saida, "resumo_metodos.csv"), row.names = FALSE)
print(resumo, digits = 4)
cat("Concluído: resultados/06_comparacao/\n")
