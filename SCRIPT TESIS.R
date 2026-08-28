# ==============================================================================
# ANALISIS ESTADISTICO ACTUALIZADO - TESIS ANDREA KARINA PEÑA
# Proyecto DIAPETICS
#
# Base integrada:
#   dataset_final_consolidado_integrado_CORREGIDO.xlsx
# Hoja:
#   dataset_final_integrado
#
# ENFOQUE ACTUAL
# --------------
# - 10 participantes, 16 sesiones, 32 mediciones pie-sesion.
# - El perfil glucemico de 7 dias se considera la exposicion principal.
# - El analisis principal usa sesiones con >=5 dias de registro glucemico.
# - La glucemia puntual de la sesion se conserva como analisis de sensibilidad.
# - La variabilidad glucemica (CV de 7 dias) se analiza de forma exploratoria.
# - PMAX bilateral es el desenlace principal.
# - PTI bilateral es un desenlace secundario.
# - ITB se mantiene descriptivo por su escasa variabilidad.
#
# IMPORTANTE
# ----------
# - Estudio piloto y exploratorio.
# - Las medidas repetidas no convierten 10 pacientes en 32 sujetos independientes.
# - No se calcula poder estadistico post hoc.
# - No se realizan afirmaciones causales ni predictivas.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACION
# ==============================================================================

ARCHIVO_EXCEL <-
  "C:/Users/Usuario/Documents/ENTREGA DE TESIS FINAL AGOSTO/Estadistica Final/dataset_final_consolidado_integrado_CORREGIDO.xlsx"
HOJA_EXCEL <- "dataset_final_integrado"

CARPETA_SALIDA <- "resultados_tesis_karina_FINAL_CORREGIDO"

DIR_TABLAS <- file.path(CARPETA_SALIDA, "tablas")
DIR_GRAFICAS <- file.path(CARPETA_SALIDA, "graficas")
DIR_MODELOS <- file.path(CARPETA_SALIDA, "modelos")

dir.create(DIR_TABLAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_GRAFICAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MODELOS, recursive = TRUE, showWarnings = FALSE)

paquetes <- c(
  "readxl", "dplyr", "tidyr", "ggplot2",
  "lme4", "lmerTest", "pbkrtest"
)

faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltantes) > 0) {
  install.packages(faltantes, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

tema_tesis <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


# ==============================================================================
# 1. FUNCIONES AUXILIARES
# ==============================================================================

exportar_csv <- function(objeto, nombre_archivo) {
  write.csv(
    objeto,
    file = file.path(DIR_TABLAS, nombre_archivo),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

resumen_numerico <- function(x, nombre_variable, unidad = "") {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(
      data.frame(
        variable = nombre_variable,
        unidad = unidad,
        n = 0,
        media = NA_real_,
        desviacion_estandar = NA_real_,
        mediana = NA_real_,
        q1 = NA_real_,
        q3 = NA_real_,
        minimo = NA_real_,
        maximo = NA_real_
      )
    )
  }
  
  data.frame(
    variable = nombre_variable,
    unidad = unidad,
    n = length(x),
    media = mean(x),
    desviacion_estandar = ifelse(length(x) > 1, sd(x), NA_real_),
    mediana = median(x),
    q1 = as.numeric(quantile(x, 0.25)),
    q3 = as.numeric(quantile(x, 0.75)),
    minimo = min(x),
    maximo = max(x),
    stringsAsFactors = FALSE
  )
}

extraer_coeficientes <- function(modelo, nombre_modelo) {
  co <- as.data.frame(summary(modelo)$coefficients)
  co$termino <- rownames(co)
  rownames(co) <- NULL
  
  nombres_esperados <- c(
    "Estimate", "Std. Error", "df",
    "t value", "Pr(>|t|)"
  )
  
  presentes <- intersect(nombres_esperados, names(co))
  
  co <- co[, c("termino", presentes), drop = FALSE]
  
  names(co) <- gsub("Estimate", "estimacion", names(co), fixed = TRUE)
  names(co) <- gsub("Std. Error", "error_estandar", names(co), fixed = TRUE)
  names(co) <- gsub("t value", "t", names(co), fixed = TRUE)
  names(co) <- gsub("Pr(>|t|)", "p", names(co), fixed = TRUE)
  
  if ("df" %in% names(co)) {
    valor_critico <- qt(0.975, df = co$df)
  } else {
    valor_critico <- rep(1.96, nrow(co))
  }
  
  co$ic95_inferior <-
    co$estimacion - valor_critico * co$error_estandar
  
  co$ic95_superior <-
    co$estimacion + valor_critico * co$error_estandar
  
  co$modelo <- nombre_modelo
  
  co
}

anova_kr_segura <- function(modelo) {
  tryCatch(
    {
      resultado <- as.data.frame(
        anova(modelo, ddf = "Kenward-Roger")
      )
      
      resultado$termino <- rownames(resultado)
      rownames(resultado) <- NULL
      resultado
    },
    error = function(e) {
      data.frame(
        termino = "No calculable",
        error = e$message
      )
    }
  )
}

metricas_modelo <- function(modelo, nombre_modelo) {
  vc <- as.data.frame(VarCorr(modelo))
  
  var_aleatoria <-
    sum(vc$vcov[vc$grp != "Residual"], na.rm = TRUE)
  
  var_residual <-
    vc$vcov[vc$grp == "Residual"][1]
  
  matriz_x <- model.matrix(modelo)
  pred_fija <- as.vector(matriz_x %*% fixef(modelo))
  var_fija <- var(pred_fija)
  
  denominador <-
    var_fija + var_aleatoria + var_residual
  
  data.frame(
    modelo = nombre_modelo,
    n_observaciones = nobs(modelo),
    n_pacientes =
      length(unique(model.frame(modelo)$paciente_id)),
    varianza_efectos_fijos = var_fija,
    varianza_aleatoria = var_aleatoria,
    varianza_residual = var_residual,
    r2_marginal = var_fija / denominador,
    r2_condicional =
      (var_fija + var_aleatoria) / denominador,
    icc =
      var_aleatoria /
      (var_aleatoria + var_residual),
    singular = isSingular(modelo, tol = 1e-4),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# f2 de Cohen para modelos lineales mixtos a partir del R2 marginal
# ------------------------------------------------------------------------------
# Se compara el modelo completo con un modelo reducido que conserva
# exactamente la misma estructura aleatoria, pero excluye el predictor fijo
# de interés. Se utiliza:
#
# f2 = (R2m_completo - R2m_reducido) / (1 - R2m_completo)
#
# Interpretacion convencional de Cohen:
# <0.02 despreciable; 0.02-<0.15 pequeno; 0.15-<0.35 moderado; >=0.35 grande.
#
# En modelos mixtos este f2 debe interpretarse como una medida exploratoria
# basada en el R2 marginal (contribucion de los efectos fijos).

calcular_f2_cohen_lmm <- function(
    modelo_completo,
    modelo_reducido,
    nombre_modelo,
    predictor
) {

  met_completo <- metricas_modelo(
    modelo_completo,
    paste0(nombre_modelo, "_completo")
  )

  met_reducido <- metricas_modelo(
    modelo_reducido,
    paste0(nombre_modelo, "_reducido")
  )

  r2_completo <- met_completo$r2_marginal[1]
  r2_reducido <- met_reducido$r2_marginal[1]

  f2 <- if (
    is.na(r2_completo) ||
    is.na(r2_reducido) ||
    r2_completo >= 1
  ) {
    NA_real_
  } else {
    (r2_completo - r2_reducido) /
      (1 - r2_completo)
  }

  interpretacion <- if (is.na(f2)) {
    NA_character_
  } else if (f2 < 0.02) {
    "Despreciable"
  } else if (f2 < 0.15) {
    "Pequeno"
  } else if (f2 < 0.35) {
    "Moderado"
  } else {
    "Grande"
  }

  data.frame(
    modelo = nombre_modelo,
    predictor = predictor,
    r2_marginal_completo = r2_completo,
    r2_marginal_reducido = r2_reducido,
    delta_r2_marginal = r2_completo - r2_reducido,
    f2_cohen = f2,
    interpretacion_f2 = interpretacion,
    stringsAsFactors = FALSE
  )
}

extraer_efecto_porcentual <- function(
    modelo,
    termino,
    nombre_modelo
) {
  
  tabla <- extraer_coeficientes(
    modelo,
    nombre_modelo
  )
  
  fila <- tabla[
    tabla$termino == termino,
    ,
    drop = FALSE
  ]
  
  if (nrow(fila) != 1) {
    return(data.frame())
  }
  
  data.frame(
    modelo = nombre_modelo,
    termino = termino,
    beta_log = fila$estimacion,
    cambio_porcentual =
      100 * (exp(fila$estimacion) - 1),
    cambio_porcentual_ic95_inferior =
      100 * (exp(fila$ic95_inferior) - 1),
    cambio_porcentual_ic95_superior =
      100 * (exp(fila$ic95_superior) - 1),
    p_satterthwaite = fila$p,
    stringsAsFactors = FALSE
  )
}

ajustar_lmm <- function(formula, datos) {
  lmer(
    formula = formula,
    data = datos,
    REML = TRUE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
}

crear_diagnosticos <- function(
    modelo,
    nombre_modelo
) {
  
  diagnostico <- data.frame(
    ajustados = fitted(modelo),
    residuos = resid(modelo),
    residuos_estandarizados =
      as.numeric(scale(resid(modelo)))
  )
  
  g_residuos <-
    ggplot(
      diagnostico,
      aes(x = ajustados, y = residuos)
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2) +
    labs(
      title =
        paste(
          "Residuos frente a valores ajustados -",
          nombre_modelo
        ),
      x = "Valores ajustados",
      y = "Residuos"
    ) +
    tema_tesis
  
  g_qq <-
    ggplot(
      diagnostico,
      aes(sample = residuos_estandarizados)
    ) +
    stat_qq(size = 2) +
    stat_qq_line() +
    labs(
      title =
        paste(
          "Grafico Q-Q de residuos -",
          nombre_modelo
        ),
      x = "Cuantiles teoricos",
      y = "Cuantiles observados"
    ) +
    tema_tesis
  
  ggsave(
    file.path(
      DIR_GRAFICAS,
      paste0(
        "diagnostico_residuos_",
        nombre_modelo,
        ".png"
      )
    ),
    g_residuos,
    width = 7,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(
      DIR_GRAFICAS,
      paste0(
        "diagnostico_qq_",
        nombre_modelo,
        ".png"
      )
    ),
    g_qq,
    width = 7,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}


# ------------------------------------------------------------------------------
# Predicciones del componente fijo para modelos con desenlace logaritmico
# ------------------------------------------------------------------------------
# Esta funcion se usa SOLO para graficar la relacion estimada por el modelo mixto.
# La linea corresponde al componente fijo (intercepto aleatorio = 0) y el IC95%
# se calcula mediante aproximacion de Wald en la escala logaritmica. Luego se
# retransforma a la escala original de PMAX o PTI.

prediccion_fija_log <- function(
    modelo,
    predictor,
    valores_predictor,
    multiplicador_x = 1
) {
  
  nuevos_datos <- data.frame(valores_predictor)
  names(nuevos_datos) <- predictor
  
  formula_fija <- reformulate(predictor)
  X <- model.matrix(formula_fija, data = nuevos_datos)
  
  beta <- fixef(modelo)[c("(Intercept)", predictor)]
  V <- as.matrix(vcov(modelo))[c("(Intercept)", predictor),
                               c("(Intercept)", predictor),
                               drop = FALSE]
  
  eta <- as.vector(X %*% beta)
  se_eta <- sqrt(diag(X %*% V %*% t(X)))
  
  data.frame(
    x_modelo = valores_predictor,
    x_original = valores_predictor * multiplicador_x,
    prediccion = exp(eta),
    ic95_inferior = exp(eta - 1.96 * se_eta),
    ic95_superior = exp(eta + 1.96 * se_eta)
  )
}


# ==============================================================================
# 2. IMPORTACION Y VALIDACION
# ==============================================================================

if (!file.exists(ARCHIVO_EXCEL)) {
  stop(
    paste0(
      "No se encontro el archivo: ",
      ARCHIVO_EXCEL
    )
  )
}

cat("\nHojas disponibles:\n")
print(excel_sheets(ARCHIVO_EXCEL))

base <- read_excel(
  ARCHIVO_EXCEL,
  sheet = HOJA_EXCEL
)

columnas_requeridas <- c(
  "paciente_id",
  "sesion_id",
  "pie",
  "sexo",
  "peso_kg",
  "talla_m",
  "imc",
  "glucemia_mgdl",
  "itb",
  "pmax_kpa",
  "pti_kpa_s",
  "region_pmax",
  "n_pasos",
  "fecha_sesion",
  "n_dias_registrados_7d",
  "n_mediciones_7d",
  "cobertura_momentos_pct",
  "glucemia_media_7d_mgdl",
  "glucemia_cv_7d_pct"
)

columnas_faltantes <-
  setdiff(
    columnas_requeridas,
    names(base)
  )

if (length(columnas_faltantes) > 0) {
  stop(
    paste(
      "Faltan columnas requeridas:",
      paste(
        columnas_faltantes,
        collapse = ", "
      )
    )
  )
}

# Las variables de momentos especificos pueden contener NA.
# Por eso no se consideran obligatorias para validar el analisis principal.
columnas_momentos <- c(
  "glucemia_antes_desayuno_7d",
  "glucemia_despues_desayuno_7d",
  "glucemia_antes_almuerzo_7d",
  "glucemia_despues_almuerzo_7d",
  "glucemia_despues_cena_7d"
)

columnas_momentos_disponibles <-
  intersect(
    columnas_momentos,
    names(base)
  )

base <- base %>%
  mutate(
    paciente_id = factor(paciente_id),
    sesion_id = as.character(sesion_id),
    
    visita =
      sub(".*_", "", sesion_id),
    
    visita =
      factor(
        visita,
        levels = c("D1", "D2")
      ),
    
    pie =
      factor(
        pie,
        levels = c(
          "Derecho",
          "Izquierdo"
        )
      ),
    
    sexo =
      factor(
        sexo,
        levels = c("F", "M")
      ),
    
    region_pmax =
      factor(region_pmax),
    
    fecha_sesion =
      as.Date(fecha_sesion),
    
    across(
      c(
        peso_kg,
        talla_m,
        imc,
        glucemia_mgdl,
        itb,
        pmax_kpa,
        pti_kpa_s,
        n_pasos,
        n_dias_registrados_7d,
        n_mediciones_7d,
        cobertura_momentos_pct,
        glucemia_media_7d_mgdl,
        glucemia_cv_7d_pct,
        all_of(columnas_momentos_disponibles)
      ),
      as.numeric
    ),
    
    imc_calculado =
      peso_kg / (talla_m^2),
    
    diferencia_imc =
      abs(imc - imc_calculado),
    
    # Glucemia puntual: cada 10 mg/dL
    glucemia_puntual_10 =
      glucemia_mgdl / 10,
    
    # Perfil glucemico medio: cada 10 mg/dL
    glucemia_media_7d_10 =
      glucemia_media_7d_mgdl / 10,
    
    # CV: cada aumento de 5 puntos porcentuales
    glucemia_cv_7d_5 =
      glucemia_cv_7d_pct / 5,
    
    perfil_glucemico_adecuado =
      if_else(
        n_dias_registrados_7d >= 5,
        "Si",
        "No"
      ),
    
    log_pmax =
      log(pmax_kpa),
    
    log_pti =
      log(pti_kpa_s)
  )

# ------------------------------------------------------------------
# Validaciones criticas
# ------------------------------------------------------------------

if (
  any(
    is.na(
      base[
        ,
        c(
          "paciente_id",
          "sesion_id",
          "pie",
          "pmax_kpa",
          "pti_kpa_s",
          "glucemia_media_7d_mgdl",
          "glucemia_cv_7d_pct"
        )
      ]
    )
  )
) {
  warning(
    paste(
      "Hay valores faltantes en variables centrales.",
      "Revise la base antes de interpretar los modelos."
    )
  )
}

if (
  any(base$pmax_kpa <= 0, na.rm = TRUE) ||
  any(base$pti_kpa_s <= 0, na.rm = TRUE)
) {
  stop(
    "PMAX y PTI deben ser mayores que cero."
  )
}

duplicados <- base %>%
  count(
    paciente_id,
    sesion_id,
    pie,
    name = "n"
  ) %>%
  filter(n > 1)

if (nrow(duplicados) > 0) {
  exportar_csv(
    duplicados,
    "control_duplicados.csv"
  )
  
  stop(
    paste(
      "Existen duplicados",
      "paciente-sesion-pie."
    )
  )
}

# Variables que deben ser constantes dentro de la sesion
inconsistencias_sesion <- base %>%
  group_by(
    paciente_id,
    sesion_id
  ) %>%
  summarise(
    n_glucemia_puntual =
      n_distinct(glucemia_mgdl),
    
    n_glucemia_media_7d =
      n_distinct(glucemia_media_7d_mgdl),
    
    n_cv_7d =
      n_distinct(glucemia_cv_7d_pct),
    
    n_dias_gluco =
      n_distinct(n_dias_registrados_7d),
    
    n_itb =
      n_distinct(itb),
    
    n_visita =
      n_distinct(visita),
    
    .groups = "drop"
  ) %>%
  filter(
    n_glucemia_puntual > 1 |
      n_glucemia_media_7d > 1 |
      n_cv_7d > 1 |
      n_dias_gluco > 1 |
      n_itb > 1 |
      n_visita > 1
  )

inconsistencias_paciente <- base %>%
  group_by(paciente_id) %>%
  summarise(
    n_sexo =
      n_distinct(sexo),
    
    n_peso =
      n_distinct(peso_kg),
    
    n_talla =
      n_distinct(talla_m),
    
    n_imc =
      n_distinct(imc),
    
    .groups = "drop"
  ) %>%
  filter(
    n_sexo > 1 |
      n_peso > 1 |
      n_talla > 1 |
      n_imc > 1
  )

exportar_csv(
  inconsistencias_sesion,
  "control_inconsistencias_sesion.csv"
)

exportar_csv(
  inconsistencias_paciente,
  "control_inconsistencias_paciente.csv"
)

if (nrow(inconsistencias_sesion) > 0) {
  stop(
    "Hay variables diferentes dentro de una misma sesion."
  )
}

if (nrow(inconsistencias_paciente) > 0) {
  stop(
    "Hay variables diferentes dentro de un mismo participante."
  )
}

control_imc <- base %>%
  distinct(
    paciente_id,
    peso_kg,
    talla_m,
    imc,
    imc_calculado,
    diferencia_imc
  ) %>%
  arrange(
    desc(diferencia_imc)
  )

exportar_csv(
  control_imc,
  "control_imc.csv"
)

# ==============================================================================
# 3. CONSTRUCCION DE UNIDADES DE ANALISIS
# ==============================================================================

# ------------------------------------------------------------------
# 3.1 Una fila por participante
# ------------------------------------------------------------------

base_participante <- base %>%
  distinct(
    paciente_id,
    sexo,
    peso_kg,
    talla_m,
    imc
  ) %>%
  arrange(paciente_id)

# ------------------------------------------------------------------
# 3.2 Una fila por paciente-sesion
# ------------------------------------------------------------------

base_sesion <- base %>%
  group_by(
    paciente_id,
    sesion_id
  ) %>%
  summarise(
    fecha_sesion =
      first(fecha_sesion),
    
    visita =
      first(visita),
    
    sexo =
      first(sexo),
    
    peso_kg =
      first(peso_kg),
    
    talla_m =
      first(talla_m),
    
    imc =
      first(imc),
    
    # Glucemia puntual previa
    glucemia_mgdl =
      first(glucemia_mgdl),
    
    glucemia_puntual_10 =
      first(glucemia_puntual_10),
    
    # Nuevo perfil glucemico
    n_dias_registrados_7d =
      first(n_dias_registrados_7d),
    
    n_mediciones_7d =
      first(n_mediciones_7d),
    
    cobertura_momentos_pct =
      first(cobertura_momentos_pct),
    
    glucemia_media_7d_mgdl =
      first(glucemia_media_7d_mgdl),
    
    glucemia_media_7d_10 =
      first(glucemia_media_7d_10),
    
    glucemia_cv_7d_pct =
      first(glucemia_cv_7d_pct),
    
    glucemia_cv_7d_5 =
      first(glucemia_cv_7d_5),
    
    perfil_glucemico_adecuado =
      first(perfil_glucemico_adecuado),
    
    # Promedios de los momentos glucemicos disponibles
    # Se conservan a nivel paciente-sesion para los descriptivos.
    glucemia_antes_desayuno_7d =
      first(glucemia_antes_desayuno_7d),
    
    glucemia_despues_desayuno_7d =
      first(glucemia_despues_desayuno_7d),
    
    glucemia_antes_almuerzo_7d =
      first(glucemia_antes_almuerzo_7d),
    
    glucemia_despues_almuerzo_7d =
      first(glucemia_despues_almuerzo_7d),
    
    glucemia_despues_cena_7d =
      first(glucemia_despues_cena_7d),
    
    # ITB
    itb =
      first(itb),
    
    # Desenlaces bilaterales
    pmax_bilateral =
      mean(pmax_kpa, na.rm = TRUE),
    
    pti_bilateral =
      mean(pti_kpa_s, na.rm = TRUE),
    
    pmax_derecho =
      pmax_kpa[pie == "Derecho"][1],
    
    pmax_izquierdo =
      pmax_kpa[pie == "Izquierdo"][1],
    
    pti_derecho =
      pti_kpa_s[pie == "Derecho"][1],
    
    pti_izquierdo =
      pti_kpa_s[pie == "Izquierdo"][1],
    
    n_pasos_derecho =
      n_pasos[pie == "Derecho"][1],
    
    n_pasos_izquierdo =
      n_pasos[pie == "Izquierdo"][1],
    
    n_pasos_min =
      min(n_pasos, na.rm = TRUE),
    
    n_pasos_total =
      sum(n_pasos, na.rm = TRUE),
    
    n_pies =
      n(),
    
    .groups = "drop"
  ) %>%
  mutate(
    log_pmax_bilateral =
      log(pmax_bilateral),
    
    log_pti_bilateral =
      log(pti_bilateral),
    
    calidad_pasos =
      if_else(
        n_pasos_min >= 3,
        "Tres o mas pasos por pie",
        "Menos de tres pasos en algun pie"
      )
  ) %>%
  arrange(
    paciente_id,
    fecha_sesion
  )

if (any(base_sesion$n_pies != 2)) {
  exportar_csv(
    base_sesion %>%
      filter(n_pies != 2),
    "control_sesiones_sin_dos_pies.csv"
  )
  
  stop(
    paste(
      "Una o mas sesiones",
      "no tienen exactamente dos pies."
    )
  )
}

# ------------------------------------------------------------------
# 3.3 Base pie-sesion
# ------------------------------------------------------------------

base_pie <- base %>%
  arrange(
    paciente_id,
    fecha_sesion,
    pie
  )

# ------------------------------------------------------------------
# 3.4 Base principal para perfiles glucemicos adecuados
#     Criterio operativo: >=5 dias registrados
# ------------------------------------------------------------------

base_sesion_principal <- base_sesion %>%
  filter(
    perfil_glucemico_adecuado == "Si"
  ) %>%
  droplevels()

base_pie_principal <- base_pie %>%
  filter(
    n_dias_registrados_7d >= 5
  ) %>%
  droplevels()

cat("\n=============================================\n")
cat("BASE PRINCIPAL DE PERFIL GLUCEMICO\n")
cat("=============================================\n")
cat(
  "Sesiones:",
  nrow(base_sesion_principal),
  "\n"
)
cat(
  "Participantes:",
  n_distinct(base_sesion_principal$paciente_id),
  "\n"
)
cat(
  "Sesiones excluidas por cobertura <5 dias:",
  nrow(base_sesion) - nrow(base_sesion_principal),
  "\n"
)

exportar_csv(
  base_participante,
  "base_participante.csv"
)

exportar_csv(
  base_sesion,
  "base_sesion_bilateral_todas.csv"
)

exportar_csv(
  base_sesion_principal,
  "base_sesion_bilateral_principal.csv"
)

exportar_csv(
  base_pie,
  "base_pie_sesion_todas.csv"
)

exportar_csv(
  base_pie_principal,
  "base_pie_sesion_principal.csv"
)


# ==============================================================================
# 4. ANALISIS DESCRIPTIVO
# ==============================================================================

# ------------------------------------------------------------------
# 4.1 Participantes
# ------------------------------------------------------------------

resumen_participantes <- bind_rows(
  resumen_numerico(
    base_participante$peso_kg,
    "Peso",
    "kg"
  ),
  
  resumen_numerico(
    base_participante$talla_m,
    "Talla",
    "m"
  ),
  
  resumen_numerico(
    base_participante$imc,
    "Indice de masa corporal",
    "kg/m2"
  )
)

sexo_participantes <- base_participante %>%
  count(sexo, name = "n") %>%
  mutate(
    porcentaje =
      100 * n / sum(n)
  )

sesiones_por_paciente <- base_sesion %>%
  count(
    paciente_id,
    name = "n_sesiones"
  ) %>%
  arrange(paciente_id)

# ------------------------------------------------------------------
# 4.2 Cobertura del perfil glucemico
# ------------------------------------------------------------------

cobertura_glucemica <- base_sesion %>%
  select(
    paciente_id,
    sesion_id,
    fecha_sesion,
    n_dias_registrados_7d,
    n_mediciones_7d,
    cobertura_momentos_pct,
    perfil_glucemico_adecuado
  ) %>%
  arrange(
    n_dias_registrados_7d,
    paciente_id
  )

resumen_cobertura <- bind_rows(
  resumen_numerico(
    base_sesion$n_dias_registrados_7d,
    "Dias registrados del perfil",
    "dias"
  ),
  
  resumen_numerico(
    base_sesion$n_mediciones_7d,
    "Mediciones glucemicas en la ventana",
    "mediciones"
  ),
  
  resumen_numerico(
    base_sesion$cobertura_momentos_pct,
    "Cobertura de momentos disponibles",
    "%"
  )
)

# ------------------------------------------------------------------
# 4.3 Variables por sesion
# ------------------------------------------------------------------

resumen_sesiones <- bind_rows(
  resumen_numerico(
    base_sesion$glucemia_mgdl,
    "Glucemia puntual de la sesion",
    "mg/dL"
  ),
  
  resumen_numerico(
    base_sesion$glucemia_media_7d_mgdl,
    "Glucemia media del perfil de 7 dias",
    "mg/dL"
  ),
  
  resumen_numerico(
    base_sesion$glucemia_cv_7d_pct,
    "Coeficiente de variacion glucemica de 7 dias",
    "%"
  ),
  
  resumen_numerico(
    base_sesion$itb,
    "Indice tobillo-brazo",
    "razon"
  ),
  
  resumen_numerico(
    base_sesion$pmax_bilateral,
    "PMAX bilateral",
    "kPa"
  ),
  
  resumen_numerico(
    base_sesion$pti_bilateral,
    "PTI bilateral",
    "kPa*s"
  ),
  
  resumen_numerico(
    base_sesion$n_pasos_total,
    "Numero total de pasos por sesion",
    "pasos"
  )
)

# Resumen de los cinco momentos disponibles, si existen
resumen_momentos <- data.frame()

if (
  length(columnas_momentos_disponibles) > 0
) {
  
  etiquetas_momentos <- c(
    glucemia_antes_desayuno_7d =
      "Antes del desayuno",
    
    glucemia_despues_desayuno_7d =
      "Despues del desayuno",
    
    glucemia_antes_almuerzo_7d =
      "Antes del almuerzo",
    
    glucemia_despues_almuerzo_7d =
      "Despues del almuerzo",
    
    glucemia_despues_cena_7d =
      "Despues de la cena"
  )
  
  resumen_momentos <- bind_rows(
    lapply(
      columnas_momentos_disponibles,
      function(v) {
        resumen_numerico(
          base_sesion[[v]],
          etiquetas_momentos[[v]],
          "mg/dL"
        )
      }
    )
  )
}

frecuencia_itb <- base_sesion %>%
  count(
    itb,
    name = "n_sesiones"
  ) %>%
  mutate(
    porcentaje =
      100 * n_sesiones /
      sum(n_sesiones)
  ) %>%
  arrange(itb)

# ------------------------------------------------------------------
# 4.4 Variables por pie
# ------------------------------------------------------------------

resumen_pie <- base_pie %>%
  group_by(pie) %>%
  summarise(
    n = n(),
    
    pmax_media =
      mean(pmax_kpa),
    
    pmax_de =
      sd(pmax_kpa),
    
    pmax_mediana =
      median(pmax_kpa),
    
    pmax_q1 =
      quantile(pmax_kpa, 0.25),
    
    pmax_q3 =
      quantile(pmax_kpa, 0.75),
    
    pti_media =
      mean(pti_kpa_s),
    
    pti_de =
      sd(pti_kpa_s),
    
    pti_mediana =
      median(pti_kpa_s),
    
    pti_q1 =
      quantile(pti_kpa_s, 0.25),
    
    pti_q3 =
      quantile(pti_kpa_s, 0.75),
    
    .groups = "drop"
  )

frecuencia_region <- base_pie %>%
  count(
    region_pmax,
    name = "n"
  ) %>%
  mutate(
    porcentaje =
      100 * n / sum(n)
  ) %>%
  arrange(desc(n))

region_por_pie <- base_pie %>%
  count(
    pie,
    region_pmax,
    name = "n"
  ) %>%
  group_by(pie) %>%
  mutate(
    porcentaje_dentro_pie =
      100 * n / sum(n)
  ) %>%
  ungroup()

# ------------------------------------------------------------------
# Exportar descriptivos
# ------------------------------------------------------------------

exportar_csv(
  resumen_participantes,
  "tabla_01_resumen_participantes.csv"
)

exportar_csv(
  sexo_participantes,
  "tabla_02_distribucion_sexo.csv"
)

exportar_csv(
  sesiones_por_paciente,
  "tabla_03_sesiones_por_paciente.csv"
)

exportar_csv(
  cobertura_glucemica,
  "tabla_04_cobertura_perfil_glucemico.csv"
)

exportar_csv(
  resumen_cobertura,
  "tabla_05_resumen_cobertura.csv"
)

exportar_csv(
  resumen_sesiones,
  "tabla_06_resumen_sesiones.csv"
)

exportar_csv(
  resumen_momentos,
  "tabla_07_resumen_momentos_glucemicos.csv"
)

exportar_csv(
  frecuencia_itb,
  "tabla_08_frecuencia_itb.csv"
)

exportar_csv(
  resumen_pie,
  "tabla_09_resumen_por_pie.csv"
)

exportar_csv(
  frecuencia_region,
  "tabla_10_region_pmax.csv"
)

exportar_csv(
  region_por_pie,
  "tabla_11_region_pmax_por_pie.csv"
)


# ==============================================================================
# 5. GRAFICAS DESCRIPTIVAS NO MODELADAS
# ==============================================================================
# Las figuras que muestran asociaciones glucemicas con PMAX/PTI se generan
# DESPUES de ajustar cada modelo, para que la linea e IC95% provengan del LMM
# y no de una regresion lineal ordinaria que ignore las medidas repetidas.

# ------------------------------------------------------------------
# 5.1 Comparacion entre pies para PMAX
# ------------------------------------------------------------------

g_pies <- base_pie %>%
  ggplot(
    aes(
      x = pie,
      y = pmax_kpa,
      group = interaction(paciente_id, sesion_id)
    )
  ) +
  geom_line(alpha = 0.45) +
  geom_point(size = 2.5) +
  labs(
    title = "PMAX del pie derecho e izquierdo por sesion",
    x = "Pie",
    y = "PMAX (kPa)"
  ) +
  tema_tesis

ggsave(
  file.path(DIR_GRAFICAS, "figura_05_pmax_por_pie.png"),
  g_pies,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------
# 5.2 Region del pico de presion
# ------------------------------------------------------------------

g_region <- ggplot(
  frecuencia_region,
  aes(
    x = reorder(region_pmax, n),
    y = n
  )
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(porcentaje, 1),
        "%)"
      )
    ),
    hjust = -0.05,
    size = 4
  ) +
  expand_limits(
    y = max(frecuencia_region$n) * 1.18
  ) +
  labs(
    title = "Distribucion de la region del pico de presion",
    x = "Region plantar",
    y = "Numero de mediciones pie-sesion"
  ) +
  tema_tesis

ggsave(
  file.path(DIR_GRAFICAS, "figura_06_region_pmax.png"),
  g_region,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 6. MODELO PRINCIPAL
#    PMAX ~ GLUCEMIA MEDIA 7 DIAS
# ==============================================================================

modelo_pmax_media7d <- ajustar_lmm(
  log_pmax_bilateral ~
    glucemia_media_7d_10 +
    (1 | paciente_id),
  datos = base_sesion_principal
)

coef_pmax_media7d <- extraer_coeficientes(
  modelo_pmax_media7d,
  "PMAX_media7d_principal"
)

metricas_pmax_media7d <- metricas_modelo(
  modelo_pmax_media7d,
  "PMAX_media7d_principal"
)

efecto_pmax_media7d <- extraer_efecto_porcentual(
  modelo_pmax_media7d,
  termino = "glucemia_media_7d_10",
  nombre_modelo = "PMAX_media7d_principal"
)

anova_pmax_media7d_kr <-
  anova_kr_segura(
    modelo_pmax_media7d
  )

exportar_csv(
  coef_pmax_media7d,
  "modelo_01_pmax_media7d_coeficientes.csv"
)

exportar_csv(
  metricas_pmax_media7d,
  "modelo_01_pmax_media7d_metricas.csv"
)

exportar_csv(
  efecto_pmax_media7d,
  "modelo_01_pmax_media7d_efecto_porcentual.csv"
)

exportar_csv(
  anova_pmax_media7d_kr,
  "modelo_01_pmax_media7d_kenward_roger.csv"
)

saveRDS(
  modelo_pmax_media7d,
  file.path(
    DIR_MODELOS,
    "modelo_pmax_media7d_principal.rds"
  )
)

crear_diagnosticos(
  modelo_pmax_media7d,
  "pmax_media7d_principal"
)

# Figura basada en el modelo mixto principal
pred_pmax_media7d <- prediccion_fija_log(
  modelo = modelo_pmax_media7d,
  predictor = "glucemia_media_7d_10",
  valores_predictor = seq(
    min(base_sesion_principal$glucemia_media_7d_10, na.rm = TRUE),
    max(base_sesion_principal$glucemia_media_7d_10, na.rm = TRUE),
    length.out = 100
  ),
  multiplicador_x = 10
)

g1 <- ggplot(
  base_sesion_principal,
  aes(
    x = glucemia_media_7d_mgdl,
    y = pmax_bilateral
  )
) +
  geom_line(
    aes(group = paciente_id),
    alpha = 0.30
  ) +
  geom_point(
    aes(shape = visita),
    size = 3
  ) +
  geom_ribbon(
    data = pred_pmax_media7d,
    aes(
      x = x_original,
      ymin = ic95_inferior,
      ymax = ic95_superior
    ),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = pred_pmax_media7d,
    aes(
      x = x_original,
      y = prediccion
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  labs(
    title = "Perfil glucemico medio de 7 dias y PMAX bilateral",
    subtitle = "Linea e IC95% derivados del modelo lineal mixto",
    x = "Glucemia media de 7 dias (mg/dL)",
    y = "PMAX bilateral (kPa)",
    shape = "Visita"
  ) +
  tema_tesis

ggsave(
  file.path(DIR_GRAFICAS, "figura_01_media7d_pmax_modelo_mixto.png"),
  g1,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 7. MODELO SECUNDARIO
#    PTI ~ GLUCEMIA MEDIA 7 DIAS
# ==============================================================================

modelo_pti_media7d <- ajustar_lmm(
  log_pti_bilateral ~
    glucemia_media_7d_10 +
    (1 | paciente_id),
  datos = base_sesion_principal
)

coef_pti_media7d <- extraer_coeficientes(
  modelo_pti_media7d,
  "PTI_media7d_secundario"
)

metricas_pti_media7d <- metricas_modelo(
  modelo_pti_media7d,
  "PTI_media7d_secundario"
)

efecto_pti_media7d <- extraer_efecto_porcentual(
  modelo_pti_media7d,
  termino = "glucemia_media_7d_10",
  nombre_modelo = "PTI_media7d_secundario"
)

anova_pti_media7d_kr <-
  anova_kr_segura(
    modelo_pti_media7d
  )

exportar_csv(
  coef_pti_media7d,
  "modelo_02_pti_media7d_coeficientes.csv"
)

exportar_csv(
  metricas_pti_media7d,
  "modelo_02_pti_media7d_metricas.csv"
)

exportar_csv(
  efecto_pti_media7d,
  "modelo_02_pti_media7d_efecto_porcentual.csv"
)

exportar_csv(
  anova_pti_media7d_kr,
  "modelo_02_pti_media7d_kenward_roger.csv"
)

saveRDS(
  modelo_pti_media7d,
  file.path(
    DIR_MODELOS,
    "modelo_pti_media7d_secundario.rds"
  )
)

crear_diagnosticos(
  modelo_pti_media7d,
  "pti_media7d_secundario"
)

# Figura basada en el modelo mixto secundario
pred_pti_media7d <- prediccion_fija_log(
  modelo = modelo_pti_media7d,
  predictor = "glucemia_media_7d_10",
  valores_predictor = seq(
    min(base_sesion_principal$glucemia_media_7d_10, na.rm = TRUE),
    max(base_sesion_principal$glucemia_media_7d_10, na.rm = TRUE),
    length.out = 100
  ),
  multiplicador_x = 10
)

g2 <- ggplot(
  base_sesion_principal,
  aes(
    x = glucemia_media_7d_mgdl,
    y = pti_bilateral
  )
) +
  geom_line(
    aes(group = paciente_id),
    alpha = 0.30
  ) +
  geom_point(
    aes(shape = visita),
    size = 3
  ) +
  geom_ribbon(
    data = pred_pti_media7d,
    aes(
      x = x_original,
      ymin = ic95_inferior,
      ymax = ic95_superior
    ),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = pred_pti_media7d,
    aes(
      x = x_original,
      y = prediccion
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  labs(
    title = "Perfil glucemico medio de 7 dias y PTI bilateral",
    subtitle = "Linea e IC95% derivados del modelo lineal mixto",
    x = "Glucemia media de 7 dias (mg/dL)",
    y = "PTI bilateral (kPa*s)",
    shape = "Visita"
  ) +
  tema_tesis

ggsave(
  file.path(DIR_GRAFICAS, "figura_02_media7d_pti_modelo_mixto.png"),
  g2,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 8. VARIABILIDAD GLUCEMICA
#    CV 7 DIAS ~ PMAX / PTI
# ==============================================================================

# ------------------------------------------------------------------
# 8.1 PMAX y CV glucemico
#     Interpretacion: por cada 5 puntos porcentuales de CV
# ------------------------------------------------------------------

modelo_pmax_cv <- ajustar_lmm(
  log_pmax_bilateral ~
    glucemia_cv_7d_5 +
    (1 | paciente_id),
  datos = base_sesion_principal
)

coef_pmax_cv <- extraer_coeficientes(
  modelo_pmax_cv,
  "PMAX_CV7d_exploratorio"
)

metricas_pmax_cv <- metricas_modelo(
  modelo_pmax_cv,
  "PMAX_CV7d_exploratorio"
)

efecto_pmax_cv <- extraer_efecto_porcentual(
  modelo_pmax_cv,
  termino = "glucemia_cv_7d_5",
  nombre_modelo =
    "PMAX_CV7d_exploratorio"
)

anova_pmax_cv_kr <-
  anova_kr_segura(
    modelo_pmax_cv
  )

exportar_csv(
  coef_pmax_cv,
  "modelo_03_pmax_cv7d_coeficientes.csv"
)

exportar_csv(
  metricas_pmax_cv,
  "modelo_03_pmax_cv7d_metricas.csv"
)

exportar_csv(
  efecto_pmax_cv,
  "modelo_03_pmax_cv7d_efecto_porcentual.csv"
)

exportar_csv(
  anova_pmax_cv_kr,
  "modelo_03_pmax_cv7d_kenward_roger.csv"
)

saveRDS(
  modelo_pmax_cv,
  file.path(
    DIR_MODELOS,
    "modelo_pmax_cv7d.rds"
  )
)

crear_diagnosticos(
  modelo_pmax_cv,
  "pmax_cv7d_exploratorio"
)

pred_pmax_cv <- prediccion_fija_log(
  modelo = modelo_pmax_cv,
  predictor = "glucemia_cv_7d_5",
  valores_predictor = seq(
    min(base_sesion_principal$glucemia_cv_7d_5, na.rm = TRUE),
    max(base_sesion_principal$glucemia_cv_7d_5, na.rm = TRUE),
    length.out = 100
  ),
  multiplicador_x = 5
)

g3 <- ggplot(
  base_sesion_principal,
  aes(
    x = glucemia_cv_7d_pct,
    y = pmax_bilateral
  )
) +
  geom_line(
    aes(group = paciente_id),
    alpha = 0.30
  ) +
  geom_point(
    aes(shape = visita),
    size = 3
  ) +
  geom_ribbon(
    data = pred_pmax_cv,
    aes(
      x = x_original,
      ymin = ic95_inferior,
      ymax = ic95_superior
    ),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = pred_pmax_cv,
    aes(
      x = x_original,
      y = prediccion
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  labs(
    title = "Variabilidad glucemica de 7 dias y PMAX bilateral",
    subtitle = "Linea e IC95% derivados del modelo lineal mixto",
    x = "CV glucemico de 7 dias (%)",
    y = "PMAX bilateral (kPa)",
    shape = "Visita"
  ) +
  tema_tesis

ggsave(
  file.path(DIR_GRAFICAS, "figura_03_cv7d_pmax_modelo_mixto.png"),
  g3,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------
# 8.2 PTI y CV glucemico
# ------------------------------------------------------------------

modelo_pti_cv <- ajustar_lmm(
  log_pti_bilateral ~
    glucemia_cv_7d_5 +
    (1 | paciente_id),
  datos = base_sesion_principal
)

coef_pti_cv <- extraer_coeficientes(
  modelo_pti_cv,
  "PTI_CV7d_exploratorio"
)

metricas_pti_cv <- metricas_modelo(
  modelo_pti_cv,
  "PTI_CV7d_exploratorio"
)

efecto_pti_cv <- extraer_efecto_porcentual(
  modelo_pti_cv,
  termino = "glucemia_cv_7d_5",
  nombre_modelo =
    "PTI_CV7d_exploratorio"
)

anova_pti_cv_kr <-
  anova_kr_segura(
    modelo_pti_cv
  )

exportar_csv(
  coef_pti_cv,
  "modelo_04_pti_cv7d_coeficientes.csv"
)

exportar_csv(
  metricas_pti_cv,
  "modelo_04_pti_cv7d_metricas.csv"
)

exportar_csv(
  efecto_pti_cv,
  "modelo_04_pti_cv7d_efecto_porcentual.csv"
)

exportar_csv(
  anova_pti_cv_kr,
  "modelo_04_pti_cv7d_kenward_roger.csv"
)

saveRDS(
  modelo_pti_cv,
  file.path(
    DIR_MODELOS,
    "modelo_pti_cv7d.rds"
  )
)

crear_diagnosticos(
  modelo_pti_cv,
  "pti_cv7d_exploratorio"
)

pred_pti_cv <- prediccion_fija_log(
  modelo = modelo_pti_cv,
  predictor = "glucemia_cv_7d_5",
  valores_predictor = seq(
    min(base_sesion_principal$glucemia_cv_7d_5, na.rm = TRUE),
    max(base_sesion_principal$glucemia_cv_7d_5, na.rm = TRUE),
    length.out = 100
  ),
  multiplicador_x = 5
)

g4 <- ggplot(
  base_sesion_principal,
  aes(
    x = glucemia_cv_7d_pct,
    y = pti_bilateral
  )
) +
  geom_line(
    aes(group = paciente_id),
    alpha = 0.30
  ) +
  geom_point(
    aes(shape = visita),
    size = 3
  ) +
  geom_ribbon(
    data = pred_pti_cv,
    aes(
      x = x_original,
      ymin = ic95_inferior,
      ymax = ic95_superior
    ),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = pred_pti_cv,
    aes(
      x = x_original,
      y = prediccion
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  labs(
    title = "Variabilidad glucemica de 7 dias y PTI bilateral",
    subtitle = "Linea e IC95% derivados del modelo lineal mixto",
    x = "CV glucemico de 7 dias (%)",
    y = "PTI bilateral (kPa*s)",
    shape = "Visita"
  ) +
  tema_tesis

ggsave(
  file.path(DIR_GRAFICAS, "figura_04_cv7d_pti_modelo_mixto.png"),
  g4,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)




# ==============================================================================
# 8.3 TAMANO DEL EFECTO - f2 DE COHEN
# ==============================================================================
# El objetivo especifico 3 solicita reportar f2 de Cohen.
# Para cada desenlace se ajusta un modelo reducido con la misma estructura
# aleatoria y sin el predictor glucemico. Los modelos completos y reducidos
# usan exactamente la misma base analitica (base_sesion_principal).

modelo_reducido_pmax <- ajustar_lmm(
  log_pmax_bilateral ~
    1 +
    (1 | paciente_id),
  datos = base_sesion_principal
)

modelo_reducido_pti <- ajustar_lmm(
  log_pti_bilateral ~
    1 +
    (1 | paciente_id),
  datos = base_sesion_principal
)

f2_pmax_media7d <- calcular_f2_cohen_lmm(
  modelo_completo = modelo_pmax_media7d,
  modelo_reducido = modelo_reducido_pmax,
  nombre_modelo = "PMAX_media7d_principal",
  predictor = "Glucemia media 7 dias (+10 mg/dL)"
)

f2_pti_media7d <- calcular_f2_cohen_lmm(
  modelo_completo = modelo_pti_media7d,
  modelo_reducido = modelo_reducido_pti,
  nombre_modelo = "PTI_media7d_secundario",
  predictor = "Glucemia media 7 dias (+10 mg/dL)"
)

f2_pmax_cv <- calcular_f2_cohen_lmm(
  modelo_completo = modelo_pmax_cv,
  modelo_reducido = modelo_reducido_pmax,
  nombre_modelo = "PMAX_CV7d_exploratorio",
  predictor = "CV glucemico 7 dias (+5 puntos porcentuales)"
)

f2_pti_cv <- calcular_f2_cohen_lmm(
  modelo_completo = modelo_pti_cv,
  modelo_reducido = modelo_reducido_pti,
  nombre_modelo = "PTI_CV7d_exploratorio",
  predictor = "CV glucemico 7 dias (+5 puntos porcentuales)"
)

tabla_f2_cohen <- bind_rows(
  f2_pmax_media7d,
  f2_pti_media7d,
  f2_pmax_cv,
  f2_pti_cv
)

exportar_csv(
  tabla_f2_cohen,
  "tabla_f2_cohen_modelos_principales.csv"
)

saveRDS(
  modelo_reducido_pmax,
  file.path(
    DIR_MODELOS,
    "modelo_reducido_pmax_intercepto_aleatorio.rds"
  )
)

saveRDS(
  modelo_reducido_pti,
  file.path(
    DIR_MODELOS,
    "modelo_reducido_pti_intercepto_aleatorio.rds"
  )
)


# ==============================================================================
# 9. ANALISIS DE SENSIBILIDAD
# ==============================================================================

# ------------------------------------------------------------------
# 9.1 Mantener ambos pies por separado
#     Predictor: glucemia media de 7 dias
# ------------------------------------------------------------------

modelo_pmax_por_pie <- ajustar_lmm(
  log_pmax ~
    glucemia_media_7d_10 +
    pie +
    (1 | paciente_id),
  datos = base_pie_principal
)

coef_pmax_por_pie <- extraer_coeficientes(
  modelo_pmax_por_pie,
  "PMAX_media7d_por_pie"
)

metricas_pmax_por_pie <- metricas_modelo(
  modelo_pmax_por_pie,
  "PMAX_media7d_por_pie"
)

efecto_pmax_por_pie <- extraer_efecto_porcentual(
  modelo_pmax_por_pie,
  termino = "glucemia_media_7d_10",
  nombre_modelo =
    "PMAX_media7d_por_pie"
)

exportar_csv(
  coef_pmax_por_pie,
  "sensibilidad_01_media7d_por_pie_coeficientes.csv"
)

exportar_csv(
  metricas_pmax_por_pie,
  "sensibilidad_01_media7d_por_pie_metricas.csv"
)

exportar_csv(
  efecto_pmax_por_pie,
  "sensibilidad_01_media7d_por_pie_efecto.csv"
)

saveRDS(
  modelo_pmax_por_pie,
  file.path(
    DIR_MODELOS,
    "modelo_pmax_media7d_por_pie.rds"
  )
)

# ------------------------------------------------------------------
# 9.2 Excluir sesiones con menos de 3 pasos en alguno de los pies
#     Se aplica sobre la base principal de cobertura glucemica adecuada.
# ------------------------------------------------------------------

base_sesion_calidad_pasos <- base_sesion_principal %>%
  filter(n_pasos_min >= 3) %>%
  droplevels()

efecto_pmax_calidad_pasos <- data.frame()
metricas_pmax_calidad_pasos <- data.frame()
anova_pmax_calidad_pasos_kr <- data.frame()

if (
  nrow(base_sesion_calidad_pasos) >= 10 &&
  n_distinct(base_sesion_calidad_pasos$paciente_id) >= 7
) {
  
  modelo_pmax_calidad_pasos <- ajustar_lmm(
    log_pmax_bilateral ~
      glucemia_media_7d_10 +
      (1 | paciente_id),
    datos = base_sesion_calidad_pasos
  )
  
  coef_pmax_calidad_pasos <- extraer_coeficientes(
    modelo_pmax_calidad_pasos,
    "PMAX_media7d_calidad_pasos"
  )
  
  metricas_pmax_calidad_pasos <- metricas_modelo(
    modelo_pmax_calidad_pasos,
    "PMAX_media7d_calidad_pasos"
  )
  
  efecto_pmax_calidad_pasos <- extraer_efecto_porcentual(
    modelo_pmax_calidad_pasos,
    termino = "glucemia_media_7d_10",
    nombre_modelo = "PMAX_media7d_calidad_pasos"
  )
  
  anova_pmax_calidad_pasos_kr <- anova_kr_segura(
    modelo_pmax_calidad_pasos
  )
  
  exportar_csv(
    coef_pmax_calidad_pasos,
    "sensibilidad_calidad_pasos_coeficientes.csv"
  )
  
  exportar_csv(
    metricas_pmax_calidad_pasos,
    "sensibilidad_calidad_pasos_metricas.csv"
  )
  
  exportar_csv(
    efecto_pmax_calidad_pasos,
    "sensibilidad_calidad_pasos_efecto.csv"
  )
  
  exportar_csv(
    anova_pmax_calidad_pasos_kr,
    "sensibilidad_calidad_pasos_kenward_roger.csv"
  )
  
  saveRDS(
    modelo_pmax_calidad_pasos,
    file.path(
      DIR_MODELOS,
      "modelo_pmax_media7d_calidad_pasos.rds"
    )
  )
  
} else {
  
  modelo_pmax_calidad_pasos <- NULL
  
  writeLines(
    "No se ajusto el modelo de calidad de pasos porque quedaron muy pocas sesiones o participantes.",
    file.path(
      DIR_MODELOS,
      "sensibilidad_calidad_pasos_no_ejecutada.txt"
    )
  )
}

# ------------------------------------------------------------------
# 9.3 Reincorporar P5 / sesion con cobertura <5 dias
#     Analisis de sensibilidad, NO principal
# ------------------------------------------------------------------

modelo_pmax_media7d_todas <- ajustar_lmm(
  log_pmax_bilateral ~
    glucemia_media_7d_10 +
    (1 | paciente_id),
  datos = base_sesion
)

coef_media7d_todas <- extraer_coeficientes(
  modelo_pmax_media7d_todas,
  "PMAX_media7d_todas_sesiones"
)

metricas_media7d_todas <- metricas_modelo(
  modelo_pmax_media7d_todas,
  "PMAX_media7d_todas_sesiones"
)

efecto_media7d_todas <- extraer_efecto_porcentual(
  modelo_pmax_media7d_todas,
  termino = "glucemia_media_7d_10",
  nombre_modelo =
    "PMAX_media7d_todas_sesiones"
)

exportar_csv(
  coef_media7d_todas,
  "sensibilidad_02_media7d_todas_coeficientes.csv"
)

exportar_csv(
  metricas_media7d_todas,
  "sensibilidad_02_media7d_todas_metricas.csv"
)

exportar_csv(
  efecto_media7d_todas,
  "sensibilidad_02_media7d_todas_efecto.csv"
)

# ------------------------------------------------------------------
# 9.4 Glucemia puntual de la sesion
#     Reproduce el enfoque del analisis previo
# ------------------------------------------------------------------

modelo_pmax_puntual <- ajustar_lmm(
  log_pmax_bilateral ~
    glucemia_puntual_10 +
    (1 | paciente_id),
  datos = base_sesion
)

coef_pmax_puntual <- extraer_coeficientes(
  modelo_pmax_puntual,
  "PMAX_glucemia_puntual"
)

metricas_pmax_puntual <- metricas_modelo(
  modelo_pmax_puntual,
  "PMAX_glucemia_puntual"
)

efecto_pmax_puntual <- extraer_efecto_porcentual(
  modelo_pmax_puntual,
  termino = "glucemia_puntual_10",
  nombre_modelo =
    "PMAX_glucemia_puntual"
)

exportar_csv(
  coef_pmax_puntual,
  "sensibilidad_03_glucemia_puntual_coeficientes.csv"
)

exportar_csv(
  metricas_pmax_puntual,
  "sensibilidad_03_glucemia_puntual_metricas.csv"
)

exportar_csv(
  efecto_pmax_puntual,
  "sensibilidad_03_glucemia_puntual_efecto.csv"
)

# ------------------------------------------------------------------
# 9.5 Ajuste exploratorio por IMC
#     Solo en perfiles con >=5 dias
# ------------------------------------------------------------------

modelo_pmax_imc <- ajustar_lmm(
  log_pmax_bilateral ~
    glucemia_media_7d_10 +
    imc +
    (1 | paciente_id),
  datos = base_sesion_principal
)

coef_pmax_imc <- extraer_coeficientes(
  modelo_pmax_imc,
  "PMAX_media7d_ajustado_IMC"
)

metricas_pmax_imc <- metricas_modelo(
  modelo_pmax_imc,
  "PMAX_media7d_ajustado_IMC"
)

efecto_pmax_imc <- extraer_efecto_porcentual(
  modelo_pmax_imc,
  termino = "glucemia_media_7d_10",
  nombre_modelo =
    "PMAX_media7d_ajustado_IMC"
)

exportar_csv(
  coef_pmax_imc,
  "sensibilidad_04_media7d_ajustado_imc_coeficientes.csv"
)

exportar_csv(
  metricas_pmax_imc,
  "sensibilidad_04_media7d_ajustado_imc_metricas.csv"
)

exportar_csv(
  efecto_pmax_imc,
  "sensibilidad_04_media7d_ajustado_imc_efecto.csv"
)

saveRDS(
  modelo_pmax_imc,
  file.path(
    DIR_MODELOS,
    "modelo_pmax_media7d_ajustado_imc.rds"
  )
)

# ------------------------------------------------------------------
# 9.6 Ajuste exploratorio por visita
# ------------------------------------------------------------------

modelo_pmax_visita <- ajustar_lmm(
  log_pmax_bilateral ~
    glucemia_media_7d_10 +
    visita +
    (1 | paciente_id),
  datos = base_sesion_principal
)

coef_pmax_visita <- extraer_coeficientes(
  modelo_pmax_visita,
  "PMAX_media7d_ajustado_visita"
)

metricas_pmax_visita <- metricas_modelo(
  modelo_pmax_visita,
  "PMAX_media7d_ajustado_visita"
)

exportar_csv(
  coef_pmax_visita,
  "sensibilidad_05_media7d_ajustado_visita_coeficientes.csv"
)

exportar_csv(
  metricas_pmax_visita,
  "sensibilidad_05_media7d_ajustado_visita_metricas.csv"
)


# ==============================================================================
# 10. LEAVE-ONE-PATIENT-OUT
#     Modelo principal con glucemia media de 7 dias
# ==============================================================================

pacientes <-
  levels(
    droplevels(
      base_sesion_principal$paciente_id
    )
  )

resultado_loo <- lapply(
  pacientes,
  function(paciente_excluido) {
    
    datos_loo <-
      base_sesion_principal %>%
      filter(
        paciente_id != paciente_excluido
      ) %>%
      droplevels()
    
    ajuste <- tryCatch(
      ajustar_lmm(
        log_pmax_bilateral ~
          glucemia_media_7d_10 +
          (1 | paciente_id),
        datos = datos_loo
      ),
      error = function(e) e
    )
    
    if (inherits(ajuste, "error")) {
      return(
        data.frame(
          paciente_excluido =
            paciente_excluido,
          
          n_sesiones =
            nrow(datos_loo),
          
          n_pacientes =
            n_distinct(
              datos_loo$paciente_id
            ),
          
          beta_log = NA_real_,
          cambio_porcentual = NA_real_,
          ic95_inferior = NA_real_,
          ic95_superior = NA_real_,
          p = NA_real_,
          singular = NA,
          error = ajuste$message,
          stringsAsFactors = FALSE
        )
      )
    }
    
    tabla <- extraer_coeficientes(
      ajuste,
      paste0(
        "sin_",
        paciente_excluido
      )
    )
    
    fila <- tabla %>%
      filter(
        termino ==
          "glucemia_media_7d_10"
      )
    
    data.frame(
      paciente_excluido =
        paciente_excluido,
      
      n_sesiones =
        nrow(datos_loo),
      
      n_pacientes =
        n_distinct(
          datos_loo$paciente_id
        ),
      
      beta_log =
        fila$estimacion,
      
      cambio_porcentual =
        100 *
        (exp(fila$estimacion) - 1),
      
      ic95_inferior =
        100 *
        (exp(fila$ic95_inferior) - 1),
      
      ic95_superior =
        100 *
        (exp(fila$ic95_superior) - 1),
      
      p =
        fila$p,
      
      singular =
        isSingular(
          ajuste,
          tol = 1e-4
        ),
      
      error = "",
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

exportar_csv(
  resultado_loo,
  "sensibilidad_06_leave_one_patient_out.csv"
)

g_loo <- ggplot(
  resultado_loo,
  aes(
    x =
      reorder(
        paciente_excluido,
        cambio_porcentual
      ),
    y = cambio_porcentual
  )
) +
  geom_hline(
    yintercept =
      efecto_pmax_media7d$cambio_porcentual,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      ymin = ic95_inferior,
      ymax = ic95_superior
    ),
    width = 0.15
  ) +
  geom_point(size = 2.5) +
  coord_flip() +
  labs(
    title =
      "Sensibilidad al excluir un participante",
    
    subtitle =
      "Cambio porcentual en PMAX por cada 10 mg/dL de glucemia media de 7 dias",
    
    x =
      "Participante excluido",
    
    y =
      "Cambio porcentual estimado"
  ) +
  tema_tesis

ggsave(
  file.path(
    DIR_GRAFICAS,
    "figura_07_leave_one_patient_out.png"
  ),
  g_loo,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ==============================================================================
# 11. ITB: DESCRIPTIVO
# ==============================================================================

n_valores_itb <-
  n_distinct(base_sesion$itb)

n_sesiones_itb_no_uno <-
  sum(
    base_sesion$itb != 1,
    na.rm = TRUE
  )

nota_itb <- paste0(
  "El ITB no se incluyo en los modelos principales. ",
  "En las ",
  nrow(base_sesion),
  " sesiones se observaron ",
  n_valores_itb,
  " valores distintos y solamente ",
  n_sesiones_itb_no_uno,
  " sesiones presentaron un ITB diferente de 1.00. ",
  "La variabilidad disponible es insuficiente para estimar ",
  "de manera estable su asociacion con PMAX o PTI."
)

writeLines(
  nota_itb,
  file.path(
    DIR_MODELOS,
    "nota_metodologica_itb.txt"
  )
)


# ==============================================================================
# 12. RESUMEN INTEGRADO DE MODELOS
# ==============================================================================

resumen_efectos <- bind_rows(
  efecto_pmax_media7d,
  efecto_pti_media7d,
  efecto_pmax_cv,
  efecto_pti_cv,
  efecto_pmax_por_pie,
  efecto_pmax_calidad_pasos,
  efecto_media7d_todas,
  efecto_pmax_puntual,
  efecto_pmax_imc
)

resumen_metricas <- bind_rows(
  metricas_pmax_media7d,
  metricas_pti_media7d,
  metricas_pmax_cv,
  metricas_pti_cv,
  metricas_pmax_por_pie,
  metricas_pmax_calidad_pasos,
  metricas_media7d_todas,
  metricas_pmax_puntual,
  metricas_pmax_imc,
  metricas_pmax_visita
)

exportar_csv(
  resumen_efectos,
  "resumen_modelos_efectos_porcentuales.csv"
)

exportar_csv(
  resumen_metricas,
  "resumen_modelos_metricas.csv"
)

exportar_csv(
  tabla_f2_cohen,
  "resumen_modelos_f2_cohen.csv"
)


# ==============================================================================
# 13. RESUMEN DE TEXTO
# ==============================================================================

archivo_resumen <-
  file.path(
    CARPETA_SALIDA,
    "RESUMEN_ANALISIS.txt"
  )

sink(archivo_resumen)

cat("============================================================\n")
cat("TESIS KARINA PENA - ANALISIS ACTUALIZADO\n")
cat("============================================================\n\n")

cat("ESTRUCTURA COMPLETA\n")
cat(
  "Participantes:",
  n_distinct(base$paciente_id),
  "\n"
)
cat(
  "Sesiones:",
  n_distinct(base$sesion_id),
  "\n"
)
cat(
  "Mediciones pie-sesion:",
  nrow(base),
  "\n\n"
)

cat("PERFIL GLUCEMICO PRINCIPAL\n")
cat(
  "Sesiones con >=5 dias:",
  nrow(base_sesion_principal),
  "\n"
)
cat(
  "Participantes en modelo principal:",
  n_distinct(
    base_sesion_principal$paciente_id
  ),
  "\n"
)
cat(
  "Sesiones con <5 dias:",
  nrow(base_sesion) -
    nrow(base_sesion_principal),
  "\n\n"
)

cat("RESUMEN DE VARIABLES POR SESION\n")
print(
  resumen_sesiones,
  row.names = FALSE
)
cat("\n")

cat("MODELO PRINCIPAL\n")
cat(
  "LOG(PMAX BILATERAL) ~ GLUCEMIA MEDIA 7 DIAS/10 + (1|PACIENTE)\n"
)
print(
  summary(modelo_pmax_media7d)
)
cat("\nKenward-Roger:\n")
print(
  anova_pmax_media7d_kr
)
cat("\nEfecto porcentual:\n")
print(
  efecto_pmax_media7d,
  row.names = FALSE
)
cat("\nMetricas:\n")
print(
  metricas_pmax_media7d,
  row.names = FALSE
)
cat("\n")

cat("MODELO SECUNDARIO PTI\n")
print(
  efecto_pti_media7d,
  row.names = FALSE
)
cat("\n")

cat("VARIABILIDAD GLUCEMICA - PMAX\n")
print(
  efecto_pmax_cv,
  row.names = FALSE
)
cat("\n")

cat("VARIABILIDAD GLUCEMICA - PTI\n")
print(
  efecto_pti_cv,
  row.names = FALSE
)
cat("\n")

cat("TAMANO DEL EFECTO - f2 DE COHEN\n")
print(
  tabla_f2_cohen,
  row.names = FALSE
)
cat("\n")

cat("SENSIBILIDAD: CALIDAD DE PASOS (>=3 PASOS POR PIE)\n")
if (nrow(efecto_pmax_calidad_pasos) > 0) {
  print(
    efecto_pmax_calidad_pasos,
    row.names = FALSE
  )
  cat("Sesiones incluidas:", nrow(base_sesion_calidad_pasos), "\n")
  cat(
    "Participantes incluidos:",
    n_distinct(base_sesion_calidad_pasos$paciente_id),
    "\n"
  )
} else {
  cat("Modelo no ejecutado por tamano insuficiente.\n")
}
cat("\n")

cat("SENSIBILIDAD: TODAS LAS SESIONES\n")
print(
  efecto_media7d_todas,
  row.names = FALSE
)
cat("\n")

cat("SENSIBILIDAD: GLUCEMIA PUNTUAL\n")
print(
  efecto_pmax_puntual,
  row.names = FALSE
)
cat("\n")

cat("SENSIBILIDAD: AJUSTE POR IMC\n")
print(
  coef_pmax_imc,
  row.names = FALSE
)
cat("\n")

cat("LEAVE-ONE-PATIENT-OUT\n")
print(
  resultado_loo,
  row.names = FALSE
)
cat("\n")

cat("NOTA ITB\n")
cat(
  nota_itb,
  "\n\n"
)

cat("ADVERTENCIAS DE INTERPRETACION\n")
cat(
  "1. El estudio es piloto y exploratorio.\n"
)
cat(
  "2. El modelo principal usa 15 sesiones con al menos 5 dias de registro glucemico.\n"
)
cat(
  "3. La glucemia media de 7 dias es un resumen del perfil reciente; no equivale a HbA1c.\n"
)
cat(
  "4. El CV glucemico se analiza de manera exploratoria.\n"
)
cat(
  "5. La glucemia puntual se conserva solo como sensibilidad.\n"
)
cat(
  "6. El ajuste por IMC debe interpretarse con cautela por el reducido numero de participantes.\n"
)
cat(
  "7. No se afirma causalidad, capacidad predictiva ni riesgo futuro de ulceracion.\n"
)
cat(
  "8. Se realizo una sensibilidad excluyendo sesiones con menos de 3 pasos en algun pie.\n"
)
cat(
  "9. No se realizo poder estadistico post hoc.\n"
)

sink()

writeLines(
  capture.output(sessionInfo()),
  file.path(
    CARPETA_SALIDA,
    "sessionInfo.txt"
  )
)


# ==============================================================================
# 14. MENSAJE FINAL
# ==============================================================================

cat("\n============================================================\n")
cat("ANALISIS ACTUALIZADO FINALIZADO\n")
cat("============================================================\n")

cat(
  "Participantes totales:",
  n_distinct(base$paciente_id),
  "\n"
)

cat(
  "Sesiones totales:",
  n_distinct(base$sesion_id),
  "\n"
)

cat(
  "Sesiones del modelo principal:",
  nrow(base_sesion_principal),
  "\n"
)

cat(
  "Participantes del modelo principal:",
  n_distinct(
    base_sesion_principal$paciente_id
  ),
  "\n"
)

cat("\nRevise primero:\n")
cat("- RESUMEN_ANALISIS.txt\n")
cat("- tablas/modelo_01_pmax_media7d_efecto_porcentual.csv\n")
cat("- tablas/modelo_03_pmax_cv7d_efecto_porcentual.csv\n")
cat("- tablas/resumen_modelos_efectos_porcentuales.csv\n")
cat("- tablas/tabla_f2_cohen_modelos_principales.csv\n")
cat("- tablas/tabla_04_cobertura_perfil_glucemico.csv\n")
cat("- tablas/tabla_07_resumen_momentos_glucemicos.csv\n")
cat("- tablas/sensibilidad_calidad_pasos_efecto.csv\n")
cat("- graficas/figura_01_media7d_pmax_modelo_mixto.png\n")
cat("- graficas/figura_03_cv7d_pmax_modelo_mixto.png\n")
cat("- graficas/figura_07_leave_one_patient_out.png\n")

cat("============================================================\n")
