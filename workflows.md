**checkM:**

```bash
checkm lineage_wf -f CheckM.txt -t 24 -x fa bins checkm_out
```

**Prokka:**

```bash
 prokka \
        --outdir prokka_out/${bin_name} \
        --prefix ${bin_name} \
        --metagenome \
        --cpus 32 \
        --force \
        ${bin}
```

**HMMER:**

```bash
for hmm in "${HMM_DIR}"/COG_*_hmm.txt; do
    hmm_base=$(basename "${hmm}" .txt)

    for faa in "${FAA_DIR}"/*.faa; do
        faa_base=$(basename "${faa}" .faa)

        tbl_out="${OUT_DIR}/tbl/${faa_base}.${hmm_base}.tbl"
        domtbl_out="${OUT_DIR}/domtbl/${faa_base}.${hmm_base}.domtbl"
        raw_out="${OUT_DIR}/raw/${faa_base}.${hmm_base}.out
      hmmsearch \
        --cpu 32 \
        --tblout "${tbl_out}" \
        --domtblout "${domtbl_out}" \
       "${hmm}" "${faa}" > "${raw_out}"
    done
done
```

**GTDB-TK**

```bash
gtdbtk classify_wf --genome_dir filtered_bins/ --out_dir gtdb-tk_out --cpus 32 -x fa
```

**Emboss:**

FORM 1

```bash
fuzzpro \
        -sequence coxl_id_hits.faa \
        -pattern "AY-x-CSFR" \
        -rformat2 excel \
        -outfile coxl1_motif.txt
```

FORM 2

```bash
fuzzpro \
        -sequence coxl_id_hits.faa \
        -pattern "AYRGAGR" \
        -rformat2 excel \
        -outfile coxl2_motif_txt
```

```bash
R-script for Wilcoxon rank-sum test

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)


getwd()
setwd(" ")

data <- read.csv(" ")   


data_long <- data %>%
  pivot_longer(
    cols      = c(CoxL, CoxM, CoxS),
    names_to  = "Gene",
    values_to = "Frequency"
  )

data_long$Gene   <- factor(data_long$Gene, levels = c("CoxL", "CoxM", "CoxS"))
data_long$Genome <- factor(data_long$Genome, levels = c("MAGs", "Reference"))


final_plot <- ggboxplot(
  data_long,
  x             = "Gene",
  y             = "Frequency",
  fill          = "Genome",
  palette       = c("gray70", "gray30"),       
  width         = 0.6,
  outlier.shape = 19,                          
  outlier.size  = 1.5,
  outlier.color = "palette",
  xlab          = "Gene",
  ylab          = "Gene Frequency (per Mb)"
) +
  stat_compare_means(
    aes(group = Genome),
    method       = "wilcox.test",
    label        = "p.format",
    label.y      = 4.2,
    size         = 4.5,
    fontface     = "bold"
  ) +
  coord_cartesian(ylim = c(0, 4.5)) +          
  stat_compare_means(
    aes(group = Genome),
    method       = "wilcox.test",
    label        = "p.format",
    label.y      = 4.2,                         
    size         = 4.5,
    fontface     = "bold"
  ) +
  labs(
    title = " "
  ) +
  theme_classic(base_size = 12) +
  theme(
    # Title
    plot.title             = element_text(face = "bold", hjust = 0.5, size = 18),
    
    # Axes
    axis.title.x           = element_text(size = 14),
    axis.title.y           = element_text(size = 14),
    axis.text.x            = element_text(size = 13, face = "bold"),
    axis.text.y            = element_text(size = 12),
    
   
    legend.position        = "right",
    legend.justification   = c(1, 1),
    legend.background      = element_rect(color = "black", fill = "white", linewidth = 0.5),
    legend.title           = element_text(face = "bold", size = 12),
    legend.text            = element_text(size = 11)
  )

print(final_plot)

ggsave("cox_boxplot.png", plot = final_plot, 
       width = 9, height = 7, dpi = 600)

R-script for Clustered column graph
library("plotly")
library("tidyr")

getwd()
setwd(" ")

data <- read.csv("---.csv")   
     final_plot <- plot_ly(
       data = data,
       x = ~Class,
       y = ~CoxL,
       type = "bar",
       name = "CoxL",
       marker = list(color = "#0072B2")
     ) %>%
     add_trace(
         y = ~CoxM,
         name = "CoxM",
         marker = list(color = "#E69F00")
       ) %>%
     add_trace(
         y = ~CoxS,
         name = "CoxS",
         marker = list(color = "#009E73")
       ) %>%
     layout(
         font = list(
             family = "Arial",       
             size = 18,              
             color = "#222222"       
           ),
         xaxis = list(
             title = list(
                 text = "<b>----</b>",  
                 font = list(size = 24, family = "Arial", color = "black")
               ),
             tickangle = -45,
             tickfont = list(size = 16, family = "Arial", color = "black"), 
             showgrid = FALSE,                         
             showline = TRUE,                          
             linecolor = "black"
           ),
         yaxis = list(
             title = list(
                 text = "<b>Cox Gene Frequency (per Mbp)</b>", 
                 font = list(size = 24, family = "Arial", color = "black")
               ),
             tickfont = list(size = 16, family = "Arial", color = "black"),
             showgrid = FALSE,                         
             showline = TRUE,                         
             linecolor = "black",
             zeroline = FALSE                          
           ),
         legend = list(
             font = list(size = 24, family = "Arial", color = "black"),
             
               bgcolor = "#FFFFFF",      
             bordercolor = "#333333",   
             borderwidth = 2,         
             orientation = "v",         
             x = 1.02,                   
             xanchor = "left",
             y = 1,
             yanchor = "top"
           ),
         
           barmode = "group",
         margin = list(t = 50, b = 90, r = 120)  
     )
 final_plot <- plotly::config(
     p = final_plot,
     toImageButtonOptions = list(
         format = "svg",
         filename = "---",
         height = 600,
         width = 900
       )
   )
 
final_plot
```