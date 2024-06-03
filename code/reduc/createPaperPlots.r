library(dplyr)
library(tidyr)
library(readr)

library(latex2exp)
library(ggsci)
library(grid)
library(gridExtra)
library(ggplot2)
library(tikzDevice)
library(stringr)
library(scales)

INCH.PER.CM <- 1 / 2.54
WIDTH <- 18.1 * INCH.PER.CM

COL.WIDTH <- 8.85 * INCH.PER.CM
BASE.SIZE <- 9
FORMAT <- "tex"

sparse <- read_csv("Erg/7m_7j_False_Stupid/ExperimentRigClasses.csv") %>%
    mutate(Method = "Sparse")
medium <- read_csv("Erg/7m_7j_False_Better/ExperimentRigClasses.csv") %>%
    mutate(Method = "Medium")
dense <- read_csv("Erg/7m_7j_False_Best/ExperimentRigClasses.csv") %>%
    mutate(Method = "Dense")

combined3 <- bind_rows(sparse, medium, dense)
combined3$ElapsedReduc <- (combined3$ElapsedTimeS - combined3$ElapsedTimeSModelGen)
write.csv(combined3, "combined.csv")

df <- read_csv("combined.csv")

OUT.PATH <- "./RPlots/"
if (!dir.exists(OUT.PATH)) {
    dir.create(OUT.PATH)
}

options(
    tikzDocumentDeclaration = "\\documentclass[conference]{IEEEtran}",
    tikzLatexPackages = c(
        getOption("tikzLatexPackages"),
        "\\usepackage{amsmath}"
    ),
    tikzSanitizeCharacters = c("%", "#"),
    tikzReplacementCharacters = c("\\%", "\\#")
)

do.save.tikz <- function(g, out.name, .width, .height) {
    tikz(str_c(OUT.PATH, out.name, ".tex"), width = .width, height = .height)
    print(g)
    dev.off()
}

do.save.tikzS <- function(g, out.name, .width, .height) {
    tikz(str_c(OUT.PATH, out.name, ".tex"),
        width = .width, height = .height,
        sanitize = TRUE
    )
    print(g)
    dev.off()
}

LFD.COLOURS <- c("black", "#E69F00", "#999999", "#009371", "#beaed4", "#ed665a", "#1f78b4", "#009371")
LFD.SHAPES <- c(15, 16, 17, 4, 5, 8, 9, 20)

theme_paper_base <- function() {
    return(theme_bw(base_size = BASE.SIZE) +
        theme(
            axis.title.x = element_text(size = BASE.SIZE),
            axis.title.y = element_text(size = BASE.SIZE),
            legend.title = element_text(size = BASE.SIZE),
            legend.position = "top",
            legend.box = "vertical",
            legend.spacing.y = unit(-0.2, "cm"),
            plot.margin = unit(c(0, 0.2, 0, 0), "cm")
        ))
}

guide_paper_base <- function() {
    return(guides(
        shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        col = guide_legend(order = 2, nrow = 1, byrow = TRUE, reverse = TRUE)
    ))
}

relable_methods <- function(orig) {
    out <- c()
    for (x in orig)
    {
        if (x == "Before") {
            out <- append(out, "None")
        } else if (x == "Sparse") {
            out <- append(out, "Sparse")
        } else if (x == "Medium") {
            out <- append(out, "Medium")
        } else if (x == "Dense") {
            out <- append(out, "Dense")
        } else {
            out <- append(out, "NA")
        }
    }
    return(out)
}

relable_strategy <- function(orig) {
    out <- c()
    for (x in orig)
    {
        if (x == "MD") {
            out <- append(out, "MD")
        } else if (x == "RMD") {
            out <- append(out, "RMD")
        } else {
            out <- append(out, "NA")
        }
    }
    return(out)
}


reduction_name <- "RMD-variant" 
point_size <- 1
al <- 0.8
x_label <- "# Variables prior to reduction"

dfh <- subset(df, Difficulty == "hard")
dfh <- dfh %>%
    filter(NoMachines != 1)


tab.vars <- dfh %>%
    select(NoMachines, NoJobs, NoVariablesBefore, NoVariablesAfter, Method) %>%
    mutate(DataType = "# Variables after reduction [log]") %>%
    rename(Data = NoVariablesAfter)

tab.runtime <- dfh %>%
    select(NoMachines, NoJobs, NoVariablesBefore, ElapsedReduc, Method) %>%
    mutate(DataType = "Runtime [s], [log]") %>%
    rename(Data = ElapsedReduc)

tab.vars.runtime <- bind_rows(tab.vars, tab.runtime)
tab.vars.runtime$Strategy <- ifelse(tab.vars.runtime$Method == "Before", "MD", "RMD")


plt.joint.vars.runtime <- ggplot(tab.vars.runtime, aes(x = NoVariablesBefore, y = Data, colour = Method, shape = Strategy)) +
    geom_point(size = point_size, alpha = al) +
    facet_grid(rows = vars(DataType), scales = "free") +
    labs(x = x_label, y = "") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "lr") +
    scale_colour_manual(reduction_name, values = LFD.COLOURS[c(-1)], labels = relable_methods) +
    scale_shape_manual(values = LFD.SHAPES[c(-1)], labels = relable_strategy) +
    guide_paper_base() +
    theme_paper_base()

ggsave(plot = plt.joint.vars.runtime, filename = "JointVarsRuntime.pdf", path = OUT.PATH, width = COL.WIDTH, height = 1.25 * COL.WIDTH)
do.save.tikzS(plt.joint.vars.runtime, "genJointVarsRuntime", COL.WIDTH, 1.25 * COL.WIDTH)


tab.reduced <- dfh %>%
    select(NoMachines, NoJobs, NoVariablesBefore, QiskitDepthr, Method, EstimatedCircDepthr) %>%
    rename(QiskitDepth = QiskitDepthr, EstimatedCircDepth = EstimatedCircDepthr)

tab.before <- dfh %>%
    select(NoMachines, NoJobs, NoVariablesBefore, QiskitDepth, EstimatedCircDepth) %>%
    mutate(Method = "Before") %>%
    distinct(.keep_all = TRUE) 
tab.qisdepth.intrgates <- bind_rows(tab.reduced, tab.before)

tab.circdepth <- tab.qisdepth.intrgates %>%
    select(NoMachines, NoJobs, NoVariablesBefore, QiskitDepth, Method) %>%
    mutate(DataType = "Circuit depth [log]") %>%
    rename(Data = QiskitDepth)

tab.intrgates <- tab.qisdepth.intrgates %>%
    select(NoMachines, NoJobs, NoVariablesBefore, EstimatedCircDepth, Method) %>%
    mutate(DataType = "# Decomposed gates [log]") %>%
    rename(Data = EstimatedCircDepth)

tab.joint.circdepth.intrGates <- rbind(tab.circdepth, tab.intrgates)
tab.joint.circdepth.intrGates$Strategy <- ifelse(tab.joint.circdepth.intrGates$Method == "Before", "MD", "RMD")

plt.joint.circdepth.intrgates <- ggplot(tab.joint.circdepth.intrGates, aes(x = NoVariablesBefore, y = Data, colour = Method, shape = Strategy)) +
    geom_point(size = point_size, alpha = al) +
    facet_grid(rows = vars(DataType), scales = "free") +
    labs(x = x_label, y = "") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "lr") +
    scale_colour_manual(reduction_name, values = LFD.COLOURS, labels = relable_methods) +
    scale_shape_manual(values = LFD.SHAPES, labels = relable_strategy) +
    guide_paper_base() +
    theme_paper_base()

ggsave(plot = plt.joint.circdepth.intrgates, filename = "JointCircDepthIntrGates.pdf", path = OUT.PATH, width = COL.WIDTH, height = 1.25 * COL.WIDTH)
do.save.tikzS(plt.joint.circdepth.intrgates, "genJointCircDepthIntrGates", COL.WIDTH, 1.25 * COL.WIDTH)


tab.reduced <- dfh %>%
    select(NoMachines, NoJobs, NoVariablesBefore, Method, density_degr1, density_degr2, density_degr3, density_degr4) %>%
    rename(density_deg1 = density_degr1, density_deg2 = density_degr2, density_deg3 = density_degr3, density_deg4 = density_degr4)

tab.before <- dfh %>%
    select(NoMachines, NoJobs, NoVariablesBefore, density_deg1, density_deg2, density_deg3, density_deg4) %>%
    mutate(Method = "Before") %>%
    distinct(.keep_all = TRUE) 
tab.densities <- bind_rows(tab.reduced, tab.before)

tab.densities <- tab.densities %>%
    pivot_longer(cols = starts_with("density"), names_to = "Type", values_to = "Density")

tab.densities$Strategy <- ifelse(tab.densities$Method == "Before", "MD", "RMD")

Type.labs <- c("$k=1$", "$k=2$", "$k=3$", "$k=4$")
names(Type.labs) <- c("density_deg1", "density_deg2", "density_deg3", "density_deg4")

plt.densities <- ggplot(tab.densities, aes(x = NoVariablesBefore, y = Density, colour = Method, shape = Strategy)) +
    geom_point(size = point_size, alpha = al) +
    facet_grid(cols = vars(Type), labeller = labeller(Type = Type.labs)) +
    labs(x = x_label, y = "Density $d_k$") +
    scale_y_continuous(labels = percent) +
    scale_colour_manual(reduction_name, values = LFD.COLOURS, labels = relable_methods) +
    scale_shape_manual(values = LFD.SHAPES, labels = relable_strategy) +
    guide_paper_base() +
    theme_paper_base() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1))

ggsave(plot = plt.densities, filename = "Densities.pdf", path = OUT.PATH, width = COL.WIDTH, height = 0.75 * COL.WIDTH)
do.save.tikzS(plt.densities, "genDensities", COL.WIDTH, 0.75 * COL.WIDTH)
