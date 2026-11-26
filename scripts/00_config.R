# =============================================================================
# 00_config.R — definições comuns aos cinco métodos (Desafio 07)
# -----------------------------------------------------------------------------
# Carregado no início de cada script com source("scripts/00_config.R").
# Executar os scripts a partir da raiz do repositório.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

# ---- Definições adotadas em TODOS os métodos --------------------------------
# Idade-índice (Ii): 6 anos. Está dentro da amplitude observada (2,6 a 6,4
# anos), próxima da idade de corte usual do eucalipto e da última medição de
# 8 das 9 parcelas, o que reduz a extrapolação na classificação.
IDADE_INDICE <- 6

# Classes de índice de local (S, em m na idade-índice): centros de 2 em 2 m
# (classes de 18-20, 20-22, 22-24 e 24-26 m), cobrindo a amplitude de S das 9
# parcelas sem desenhar curvas muito fora dos dados.
# Cada classe vai de S - 1 a S + 1 m.
CLASSES_S <- c(19, 21, 23, 25)

# Grade de idades para desenhar o feixe de curvas (amplitude observada)
IDADES_CURVA <- seq(2.5, 6.5, by = 0.1)

# ---- Leitura e padronização dos dados ----------------------------------------
le_dados <- function(arquivo = "dados/Sitio.xls") {
  d <- as.data.frame(read_excel(arquivo, sheet = "Dados"))
  names(d) <- c("projeto", "talhao", "parcela", "idade", "area_basal", "hd", "volume")
  d$id <- sprintf("T%d-P%d", d$talhao, d$parcela)
  d <- d[order(d$id, d$idade), ]
  rownames(d) <- NULL
  d
}

# ---- Modelo de crescimento de base: Chapman-Richards --------------------------
#   Hd = b0 * (1 - exp(-b1 * I))^b2
chapman <- function(I, b0, b1, b2) b0 * (1 - exp(-b1 * I))^b2

# ---- Métricas de ajuste (na escala de Hd, em m) -------------------------------
metricas <- function(obs, est, p, nome = "") {
  res <- obs - est
  n <- length(obs)
  sqr <- sum(res^2)
  data.frame(
    modelo  = nome,
    n       = n,
    p       = p,
    r2      = 1 - sqr / sum((obs - mean(obs))^2),
    syx_m   = sqrt(sqr / (n - p)),
    syx_pct = 100 * sqrt(sqr / (n - p)) / mean(obs),
    vies_m  = mean(res),
    aic     = n * log(sqr / n) + 2 * (p + 1)
  )
}

# ---- Classe de S (centro mais próximo) ----------------------------------------
classe_s <- function(s) {
  sapply(s, function(x) CLASSES_S[which.min(abs(CLASSES_S - x))])
}

# ---- Medição de cada parcela mais próxima da idade-índice ---------------------
# Critério de classificação usado em todos os métodos: projeta-se para a
# idade-índice a medição de cada parcela mais próxima de Ii, que é a que exige
# menor extrapolação ao longo da curva.
medicao_mais_proxima <- function(d, ii = IDADE_INDICE) {
  do.call(rbind, lapply(split(d, d$id), function(x) x[which.min(abs(x$idade - ii)), ]))
}

# ---- Tema e cores dos gráficos ------------------------------------------------
tema <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "right",
        plot.title = element_text(face = "bold"))

# Gráfico padrão: feixe de curvas sobre as trajetórias observadas
grafico_feixe <- function(curvas, d, titulo, subtitulo = NULL) {
  ggplot() +
    geom_line(data = curvas, aes(idade, hd, group = factor(S)),
              colour = "grey35", linewidth = 0.8) +
    geom_text(data = curvas[curvas$idade == max(curvas$idade), ],
              aes(idade + 0.05, hd, label = paste0("S = ", S)),
              hjust = 0, size = 3.3, colour = "grey25") +
    geom_line(data = d, aes(idade, hd, colour = id), alpha = 0.7) +
    geom_point(data = d, aes(idade, hd, colour = id), size = 2) +
    geom_vline(xintercept = IDADE_INDICE, linetype = "dashed", colour = "firebrick") +
    scale_x_continuous(limits = c(2.4, 7.1), breaks = 2:7) +
    labs(x = "Idade (anos)", y = "Altura dominante (m)", colour = "Parcela",
         title = titulo,
         subtitle = if (is.null(subtitulo)) sprintf("Idade-índice = %g anos (linha tracejada)", IDADE_INDICE) else subtitulo) +
    tema
}

grafico_residuos <- function(df, titulo) {
  # df: colunas modelo, estimado, residuo
  ggplot(df, aes(estimado, residuo)) +
    geom_hline(yintercept = 0, colour = "firebrick") +
    geom_point(colour = "#2E7D32", alpha = 0.8) +
    facet_wrap(~ modelo) +
    labs(x = "Hd estimada (m)", y = "Resíduo (m)", title = titulo) +
    tema
}

pasta_saida <- function(nome) {
  p <- file.path("resultados", nome)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}
