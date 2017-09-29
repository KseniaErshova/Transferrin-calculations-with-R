library(ggplot2)
library(magrittr)
library(tidyr)
library(devtools)
library(gridExtra)
library(dplyr)
library(reshape2)
library (plyr)

# Stat smooth function
stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                             list(a = format(coef(m)[1], digits = 3), 
                                                  b = format(coef(m)[2], digits = 3), 
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
)



# Export data for Standard curve

setwd("~/path")      #CHANGE HERE
path = "~/path"   #CHANGE HERE
out.file<-""
file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i], skip = 8)
  out.file <- rbind(out.file, file)
}
out.file <- as.data.frame(sapply(out.file, as.numeric))
out.file = out.file[-1,]
rownames(out.file) = NULL


OD = c(0.005,	0.0025,	0.00125,	0.000625,	0.0003125,
       0.00015625,	0.000078125,	0.0000390625,	0.00001953125,	0.000009765625,	0.0000048828125, 0)
colnames(out.file) = OD

d <- melt(out.file)
d = d[rev(rownames(d)),]
colnames(d) = c("OD", "Fluorescence")
rownames(d) = NULL

# Find and subtract mean fluorescence of RPMI (baseline)
A = mean(d[1:8,2])
d$Fluorescence <- (d$Fluorescence - A)

d = d[9:length(d$Fluorescence),]
d$OD = as.numeric(as.character(d$OD))

# Mathematical transformation
d2 = d
d2$OD <- d2$OD^(0.125)

S = ggplot(d2, aes(x = OD, y = Fluorescence)) +
        geom_point() +
        xlab(expression(sqrt(OD, 8))) +
        scale_y_log10(name="Fluorescence (AU)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(250, 250000)) +
        annotation_logticks(sides = "l") +
        stat_smooth_func(data = d2, geom="text",method="lm",hjust=0,parse=T) + 
        geom_smooth(data = d2, method = "lm", se = F, colour = "black") +
        theme(panel.grid.major = element_line(linetype = "dotted", colour="grey26"), 
                panel.grid.major.x = element_blank(),
                axis.title.x = element_text(colour="grey26",size=16),
                axis.title.y = element_text(colour="grey26",size=14),
                plot.title = element_text(colour="grey26", size = 20, face = "bold.italic"),
                axis.line.x = element_line(colour = "grey26", size = 0.6),
                axis.line.y = element_line(colour = "grey26", size = 0.6),
                panel.background = element_blank()) +
        annotate("text", x=0.4, y=1000, size=8, colour="black",
                 label= "OD==(frac(log[10](Fluorescence)-2.49, 6.25))^{8}", parse = TRUE) +
        ggtitle('Standard curve Fluorescence vs OD; A.baumannii AB074')
S
ggsave("Microb_Standard_Curve.pdf", width=15, height=7, dpi=200)

# Function to test equation
p <- function(x) { ((log(x, 10) - 2.49) / 6.25)^8 }

# Add descriptive statistics
DS = summary(out.file)
capture.output(DS, file="Microb.txt", append=TRUE)

pvalues = pairwise.t.test(d$Fluorescence, d$OD, p.adjust = "none")
capture.output(pvalues, file="Microb.txt", append=TRUE)


