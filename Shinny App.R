# ============================================================
# OPV Causal ML - Interactive Design Engine
# Shiny R Application
# Author: Samuel Bizimana | JUNIA ISEN
# Supervisor: Dr. Kekeli N'KONOU
# Run with: shiny::runApp("shiny_app.R")
# ============================================================

library(shiny)
library(ggplot2)
library(data.table)
library(DT)

# ── Load data ─────────────────────────────────────────────
# Set your project path here
BASE <- "C:/Users/Samuel Bizimana/OneDrive/Desktop/Research Training"

load_csv <- function(rel_path) {
  path <- file.path(BASE, rel_path)
  if (file.exists(path)) fread(path) else NULL
}

dml_results  <- load_csv("results/tables/dml_ewg_corrected.csv")
val_results  <- load_csv("results/tables/phase3_xtb_validation_final.csv")
cates        <- load_csv("data/processed/experimental_with_cates.csv")
scaffold_val <- load_csv("results/tables/validation_by_scaffold.csv")
sens_results <- load_csv("results/tables/sensitivity_results.csv")
boot_results <- load_csv("results/tables/bootstrap_leakage_gap.csv")
nuisance     <- load_csv("results/tables/nuisance_performance.csv")

# ── Figure paths ──────────────────────────────────────────
fig <- function(name) file.path(BASE, "results/figures", name)

# ── Hammett sigma-para values ─────────────────────────────
HAMMETT <- data.frame(
  group = c("Nitro", "Cyano", "Sulfonyl", "Carbonyl", "Ester", "Halogen"),
  sigma = c(0.78, 0.66, 0.72, 0.50, 0.45, 0.23),
  stringsAsFactors = FALSE
)

# ── UI ────────────────────────────────────────────────────
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', sans-serif; background: #f8f9fa; }
      .navbar { background: #1F3A5F !important; }
      .navbar-brand, .navbar-nav > li > a {
        color: #ffffff !important; font-weight: 500;
      }
      .nav-tabs > li.active > a {
        color: #1F3A5F !important; font-weight: 600;
        border-top: 3px solid #1F3A5F;
      }
      .well { background: #ffffff; border: 1px solid #dee2e6;
              border-radius: 8px; padding: 18px; }
      .result-box {
        background: #ffffff; border: 1px solid #dee2e6;
        border-radius: 8px; padding: 16px; margin-bottom: 14px;
      }
      .result-box h5 { color: #1F3A5F; font-weight: 600;
                       margin-bottom: 8px; font-size: 13px;
                       text-transform: uppercase; letter-spacing: 0.05em; }
      .big-num { font-size: 26px; font-weight: 700;
                 color: #1F3A5F; font-family: 'Courier New', monospace; }
      .unit { font-size: 12px; color: #6c757d; margin-top: 2px; }
      .pass-badge { display: inline-block; padding: 3px 10px;
                    border-radius: 12px; font-size: 12px; font-weight: 500; }
      .pass { background: #d4edda; color: #155724; }
      .fail { background: #f8d7da; color: #721c24; }
      .warn { background: #fff3cd; color: #856404; }
      .info-box { background: #e8f0fe; border-left: 4px solid #1F3A5F;
                  padding: 12px 16px; border-radius: 0 8px 8px 0;
                  margin-bottom: 12px; font-size: 13px; }
      .section-head { font-size: 15px; font-weight: 700;
                      color: #1F3A5F; border-bottom: 2px solid #1F3A5F;
                      padding-bottom: 6px; margin-bottom: 16px; }
      .step-box { background: #f0f4ff; border: 1px solid #c5d0f0;
                  border-radius: 8px; padding: 14px; margin-bottom: 10px; }
      .step-num { font-size: 11px; font-weight: 700; color: #1F3A5F;
                  text-transform: uppercase; letter-spacing: 0.08em; }
      .mono { font-family: 'Courier New', monospace; }
      hr { border-color: #dee2e6; }
    "))
  ),
  
  navbarPage(
    title = "OPV Causal ML - Design Engine",
    id    = "nav",
    
    # ── TAB 1: Overview ──────────────────────────────────
    tabPanel("Overview",
             br(),
             fluidRow(
               column(12,
                      div(class = "well",
                          h3("Causal Machine Learning for Organic Acceptor Molecular Design",
                             style = "color:#1F3A5F; font-weight:700;"),
                          p("Samuel Bizimana | JUNIA ISEN | Supervisor: Dr. Kekeli N'KONOU",
                            style = "color:#6c757d; margin-bottom:16px;"),
                          p("This application demonstrates the causal ML pipeline applied to
               15,529 organic acceptor molecules. Select a molecule on the
               Design Engine tab to see per-molecule CATE estimates and
               counterfactual design reasoning.")
                      )
               )
             ),
             fluidRow(
               column(3,
                      div(class = "result-box",
                          h5("Primary Finding"),
                          div(class = "big-num", "-0.025"),
                          div(class = "unit", "eV - EWG effect on LUMO"),
                          br(),
                          span(class = "pass-badge pass", "p = 0.0008")
                      )
               ),
               column(3,
                      div(class = "result-box",
                          h5("HOMO Effect"),
                          div(class = "big-num", "-0.012"),
                          div(class = "unit", "eV - EWG effect on HOMO"),
                          br(),
                          span(class = "pass-badge fail", "p = 0.115 NS")
                      )
               ),
               column(3,
                      div(class = "result-box",
                          h5("xTB Validation MAE"),
                          div(class = "big-num", "0.2155"),
                          div(class = "unit", "eV - below 0.25 eV threshold"),
                          br(),
                          span(class = "pass-badge pass", "6/10 pass")
                      )
               ),
               column(3,
                      div(class = "result-box",
                          h5("Bootstrap Leakage"),
                          div(class = "big-num", "-0.002"),
                          div(class = "unit", "Mean R² gap - 20 seeds"),
                          br(),
                          span(class = "pass-badge pass", "No leakage")
                      )
               )
             ),
             fluidRow(
               column(6,
                      div(class = "well",
                          div(class = "section-head", "DML Causal Estimates"),
                          tableOutput("dml_table")
                      )
               ),
               column(6,
                      div(class = "well",
                          div(class = "section-head", "Nuisance Model Diagnostics"),
                          tableOutput("nuisance_table"),
                          br(),
                          div(class = "info-box",
                              "R²(ml_l) = 0.181 means confounders explain only 18% of HOMO
               variance. DML estimates are valid but carry higher uncertainty
               than ideal. HOMO is governed by factors not captured by
               single-molecule structural descriptors."
                          )
                      )
               )
             ),
             fluidRow(
               column(12,
                      div(class = "well",
                          div(class = "section-head", "Study Pipeline"),
                          fluidRow(
                            column(2, div(class = "step-box", div(class="step-num","Phase 0"),
                                          strong("Environment"), br(),
                                          tags$small("R 4.6.0 + Python 3.10"))),
                            column(2, div(class = "step-box", div(class="step-num","Phase 1"),
                                          strong("Dataset"), br(),
                                          tags$small("15,529 molecules"))),
                            column(2, div(class = "step-box", div(class="step-num","Phase 2"),
                                          strong("Causal Analysis"), br(),
                                          tags$small("DML + Forests + Sensitivity"))),
                            column(2, div(class = "step-box", div(class="step-num","Phase 3"),
                                          strong("Design Engine"), br(),
                                          tags$small("CATE-driven SMARTS"))),
                            column(2, div(class = "step-box", div(class="step-num","Phase 4"),
                                          strong("xTB Validation"), br(),
                                          tags$small("MAE = 0.2155 eV"))),
                            column(2, div(class = "step-box", div(class="step-num","Output"),
                                          strong("Manuscript"), br(),
                                          tags$small("Nature Machine Intelligence")))
                          )
                      )
               )
             )
    ),
    
    # ── TAB 2: Design Engine ─────────────────────────────
    tabPanel("Design Engine",
             br(),
             sidebarLayout(
               sidebarPanel(
                 div(class = "section-head", "Step 1: Select molecule"),
                 p("Choose from 598 experimental organic acceptors with
             pre-computed per-molecule CATE estimates from the
             causal forest.", style = "font-size:13px; color:#6c757d;"),
                 selectInput("mol_id", "Molecule",
                             choices  = if (!is.null(cates))
                               setNames(seq_len(nrow(cates)),
                                        paste0("Mol ", seq_len(nrow(cates)),
                                               " | HOMO: ",
                                               round(cates$homo_ev, 3),
                                               " eV | CATE: ",
                                               round(cates$CATE, 4))) else 1:10,
                             selected = 1
                 ),
                 hr(),
                 div(class = "section-head", "Step 2: Set target"),
                 sliderInput("target_shift",
                             "Target HOMO shift (eV)",
                             min = -0.5, max = 0.1, value = -0.2, step = 0.01
                 ),
                 hr(),
                 actionButton("run_btn", "Run Design Engine",
                              class = "btn btn-primary btn-block",
                              style = "background:#1F3A5F; border-color:#1F3A5F;
                     font-weight:600; width:100%;"
                 ),
                 br(),
                 div(class = "info-box",
                     strong("Population-level DML estimates:"), br(),
                     "EWG to LUMO: -0.025 eV (p=0.0008)", br(),
                     "EWG to HOMO: -0.012 eV (p=0.115)", br(),
                     "Bandgap: -0.018 eV (p=0.061)"
                 )
               ),
               
               mainPanel(
                 conditionalPanel("input.run_btn == 0",
                                  div(class = "info-box",
                                      style = "margin-top:40px; text-align:center;",
                                      h4("Select a molecule and click Run Design Engine",
                                         style = "color:#1F3A5F;"),
                                      p("The engine uses per-molecule CATE from the causal
                 forest to compute required EWG additions for your
                 target HOMO shift.")
                                  )
                 ),
                 conditionalPanel("input.run_btn > 0",
                                  uiOutput("design_output")
                 )
               )
             )
    ),
    
    # ── TAB 3: Causal Analysis ───────────────────────────
    tabPanel("Causal Analysis",
             br(),
             fluidRow(
               column(6,
                      div(class = "well",
                          div(class = "section-head", "DML Average Treatment Effects"),
                          plotOutput("dml_plot", height = "300px"),
                          br(),
                          div(class = "info-box",
                              "LUMO effect survives Bonferroni correction (p=0.0024).
               HOMO and Bandgap are not significant after correction."
                          )
                      )
               ),
               column(6,
                      div(class = "well",
                          div(class = "section-head", "Sensitivity Analysis"),
                          p("Robustness values for LUMO (primary) and HOMO (secondary):",
                            style = "font-size:13px;"),
                          fluidRow(
                            column(6,
                                   div(class = "result-box",
                                       h5("LUMO - Primary"),
                                       div(class = "big-num", "0.124"),
                                       div(class = "unit",
                                           "Robustness value - moderately sensitive"),
                                       br(),
                                       span(class = "pass-badge warn",
                                            "12.4% needed to nullify")
                                   )
                            ),
                            column(6,
                                   div(class = "result-box",
                                       h5("HOMO - Secondary"),
                                       div(class = "big-num", "0.146"),
                                       div(class = "unit", "Robustness value"),
                                       br(),
                                       span(class = "pass-badge warn",
                                            "14.6% needed to nullify")
                                   )
                            )
                          ),
                          div(class = "info-box",
                              "Both values indicate moderate sensitivity to unmeasured
               confounding. The LUMO finding should be interpreted with
               appropriate caution. Benchmark covariate: molecular weight."
                          )
                      )
               )
             ),
             fluidRow(
               column(12,
                      div(class = "well",
                          div(class = "section-head", "CATE Distribution"),
                          fluidRow(
                            column(4,
                                   div(class = "result-box",
                                       h5("Mean CATE"),
                                       div(class = "big-num", "-0.022"),
                                       div(class = "unit", "eV per unit EWG score")
                                   )
                            ),
                            column(4,
                                   div(class = "result-box",
                                       h5("CATE Range"),
                                       div(class = "big-num", "0.073"),
                                       div(class = "unit", "eV (from -0.046 to +0.027)")
                                   )
                            ),
                            column(4,
                                   div(class = "result-box",
                                       h5("Heterogeneity Test"),
                                       div(class = "big-num", "0.295"),
                                       div(class = "unit", "p-value - not confirmed"),
                                       br(),
                                       span(class = "pass-badge fail",
                                            "Not statistically significant")
                                   )
                            )
                          ),
                          div(class = "info-box",
                              "No statistically significant heterogeneity was detected
               (p=0.295). This should not be interpreted as confirming
               the absence of heterogeneity - the test may be underpowered
               at n=598. The apparent CATE-molecular weight trend in
               Figure 4 should be treated as exploratory only."
                          )
                      )
               )
             )
    ),
    
    # ── TAB 4: IV Analysis ───────────────────────────────
    tabPanel("IV Analysis",
             br(),
             div(class = "well",
                 div(class = "section-head",
                     "Instrumental Variable Analysis - No Valid Instrument Found"),
                 div(class = "info-box",
                     strong("Conclusion:"),
                     " No valid instrument was identified. Both candidates
            violate the exclusion restriction. The DML estimates
            remain the primary causal evidence."
                 ),
                 br(),
                 fluidRow(
                   column(6,
                          div(class = "step-box",
                              div(class = "step-num", "Candidate 1 - REJECTED"),
                              strong("Halogen count"), br(), br(),
                              p("Halogen count is algebraically a component of the
                 Hammett-weighted EWG score (contributing halogen x 0.23).
                 It cannot satisfy the exclusion restriction independently
                 of the treatment definition.", style="font-size:13px;")
                          )
                   ),
                   column(6,
                          div(class = "step-box",
                              div(class = "step-num", "Candidate 2 - REJECTED"),
                              strong("Non-aromatic C=C bond count"), br(), br(),
                              p("Non-aromatic C=C bonds directly affect conjugation length,
                 which is a core determinant of frontier orbital energies.
                 The instrument affects the outcome through a pathway other
                 than the EWG treatment.", style="font-size:13px;")
                          )
                   )
                 ),
                 br(),
                 fluidRow(
                   column(4,
                          div(class = "result-box",
                              h5("F-statistic"),
                              div(class = "big-num", "216.6"),
                              div(class = "unit", "Strong instrument strength")
                          )
                   ),
                   column(4,
                          div(class = "result-box",
                              h5("Wu-Hausman p-value"),
                              div(class = "big-num", "0.433"),
                              div(class = "unit",
                                  "Reported descriptively only - instrument invalid")
                          )
                   ),
                   column(4,
                          div(class = "result-box",
                              h5("Sargan test"),
                              div(class = "big-num", "N/A"),
                              div(class = "unit",
                                  "Just-identified - not applicable")
                          )
                   )
                 ),
                 div(class = "info-box",
                     "The Wu-Hausman test returned p=0.433. However, because
           neither proposed instrument could be defended under the
           exclusion restriction, this test was not treated as
           independent evidence for the consistency of OLS."
                 )
             )
    ),
    
    # ── TAB 5: xTB Validation ────────────────────────────
    tabPanel("xTB Validation",
             br(),
             fluidRow(
               column(4,
                      div(class = "result-box",
                          h5("Overall MAE"),
                          div(class = "big-num", "0.2155"),
                          div(class = "unit", "eV - PASSED 0.25 eV threshold")
                      )
               ),
               column(4,
                      div(class = "result-box",
                          h5("Candidates Passing"),
                          div(class = "big-num", "6/10"),
                          div(class = "unit", "Within 0.25 eV of xTB")
                      )
               ),
               column(4,
                      div(class = "result-box",
                          h5("All Synthesizable"),
                          div(class = "big-num", "10/10"),
                          div(class = "unit", "SAScore below 4.0")
                      )
               )
             ),
             fluidRow(
               column(7,
                      div(class = "well",
                          div(class = "section-head",
                              "ML Predicted vs xTB Computed HOMO"),
                          plotOutput("xtb_scatter", height = "350px")
                      )
               ),
               column(5,
                      div(class = "well",
                          div(class = "section-head", "Scaffold-Level Results"),
                          tableOutput("scaffold_table"),
                          br(),
                          div(class = "info-box",
                              "All three scaffolds produced at least one passing candidate.
               Scaffold HOMO=-5.99 eV shows best agreement (mean error
               0.076 eV). Scaffold HOMO=-5.65 eV shows largest discrepancy
               (mean error 0.325 eV) - CATE underestimates quantum-chemical
               response for this class."
                          )
                      )
               )
             ),
             fluidRow(
               column(12,
                      div(class = "well",
                          div(class = "section-head", "Full Validation Results"),
                          DTOutput("xtb_table")
                      )
               )
             )
    ),
    
    # ── TAB 6: About ─────────────────────────────────────
    tabPanel("About",
             br(),
             div(class = "well",
                 div(class = "section-head", "About This Project"),
                 fluidRow(
                   column(6,
                          h5("Research details"),
                          tags$ul(
                            tags$li("Learner: Samuel Bizimana"),
                            tags$li("Supervisor: Dr. Kekeli N'KONOU"),
                            tags$li("Institution: JUNIA ISEN"),
                            tags$li("Target journal: Nature Machine Intelligence"),
                            tags$li(tags$a(
                              "GitHub repository",
                              href = "https://github.com/Bsamuel-tech/opv-causal-ml",
                              target = "_blank"
                            ))
                          ),
                          br(),
                          h5("Methods used"),
                          tags$ul(
                            tags$li("Double Machine Learning (DoubleML R package)"),
                            tags$li("Causal forests (grf R package, 4,000 trees)"),
                            tags$li("Sensemakr sensitivity analysis"),
                            tags$li("RDKit molecular processing"),
                            tags$li("xTB 6.7.1 GFN2 quantum chemistry")
                          )
                   ),
                   column(6,
                          h5("Key results"),
                          tags$ul(
                            tags$li("EWG causally lowers LUMO by 0.025 eV (p=0.0008)"),
                            tags$li("HOMO not significant after confounder control (p=0.115)"),
                            tags$li("LUMO survives Bonferroni correction (p=0.0024)"),
                            tags$li("No valid IV instrument identified"),
                            tags$li("Heterogeneity not confirmed (p=0.295)"),
                            tags$li("xTB MAE=0.2155 eV, 6/10 candidates validated"),
                            tags$li("Bootstrap leakage gap mean=-0.0017 (20 seeds)")
                          ),
                          br(),
                          h5("Hammett sigma-para values"),
                          tableOutput("hammett_table")
                   )
                 )
             )
    )
  )
)

# ── Server ────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # DML table
  output$dml_table <- renderTable({
    data.frame(
      Outcome   = c("HOMO", "LUMO", "Bandgap"),
      `ATE (eV)`= c(-0.012, -0.025, -0.018),
      `P-value` = c(0.115, 0.0008, 0.061),
      `P (Bonferroni)` = c(0.345, 0.0024, 0.182),
      Significant = c("No", "Yes", "No"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Nuisance table
  output$nuisance_table <- renderTable({
    data.frame(
      Model = c("ml_l: HOMO ~ confounders",
                "ml_m: EWG ~ confounders"),
      `R²`  = c(0.181, 0.713),
      RMSE  = c("0.133 eV", "0.712 units"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # DML plot
  output$dml_plot <- renderPlot({
    df_plot <- data.frame(
      outcome = c("HOMO", "LUMO", "Bandgap"),
      ATE     = c(-0.012, -0.025, -0.018),
      SE      = c(0.0076, 0.0074, 0.0095),
      sig     = c("Not significant", "Significant", "Not significant")
    )
    df_plot$lower <- df_plot$ATE - 1.96 * df_plot$SE
    df_plot$upper <- df_plot$ATE + 1.96 * df_plot$SE
    df_plot$color <- c("#748CAB", "#2166ac", "#748CAB")
    
    ggplot(df_plot, aes(x = outcome, y = ATE, color = sig)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = "gray60") +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    width = 0.15, linewidth = 1.2) +
      geom_point(size = 5) +
      scale_color_manual(
        values = c("Significant" = "#2166ac",
                   "Not significant" = "#748CAB"),
        name = NULL
      ) +
      labs(x = NULL,
           y = "Average Treatment Effect (eV)",
           caption = "Error bars = 95% CI") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom",
            plot.caption = element_text(color = "gray50"))
  })
  
  # xTB scatter
  output$xtb_scatter <- renderPlot({
    val <- data.frame(
      candidate = 1:10,
      ml  = c(-5.9719,-5.9719,-5.9719,
              -5.6331,-5.6331,-5.6331,
              -6.1831,-5.6331,-6.1831,-5.6331),
      xtb = c(-5.8796,-5.9140,-5.8949,
              -6.1266,-6.0983,-6.1346,
              -5.9269,-5.6874,-6.1366,-5.7435),
      pass = c(TRUE,TRUE,TRUE,
               FALSE,FALSE,FALSE,
               FALSE,TRUE,TRUE,TRUE)
    )
    val$label  <- paste0("C", val$candidate)
    val$result <- ifelse(val$pass, "Pass", "Fail")
    rng <- range(c(val$ml, val$xtb)) + c(-0.05, 0.05)
    
    ggplot(val, aes(x = ml, y = xtb,
                    color = result, label = label)) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray50") +
      geom_point(size = 5, alpha = 0.9) +
      geom_text(nudge_x = 0.015, nudge_y = 0.015,
                size = 3.5, color = "gray30") +
      scale_color_manual(
        values = c("Pass" = "#2E75B6", "Fail" = "#748CAB"),
        name   = "Validation"
      ) +
      annotate("text",
               x = rng[1] + 0.05,
               y = rng[2] - 0.03,
               label = "MAE = 0.2155 eV\n6/10 pass",
               size = 4, hjust = 0,
               color = "gray20", fontface = "bold") +
      coord_cartesian(xlim = rng, ylim = rng) +
      labs(x = "ML Predicted HOMO (eV)",
           y = "xTB Calibrated HOMO (eV)",
           caption = "Dashed = perfect agreement") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
  })
  
  # Scaffold table
  output$scaffold_table <- renderTable({
    data.frame(
      Scaffold = c("HOMO=-5.99", "HOMO=-5.65", "HOMO=-6.20"),
      Candidates = c(3, 5, 2),
      Pass = c(3, 2, 1),
      Fail = c(0, 3, 1),
      `Mean error (eV)` = c(0.076, 0.325, 0.151),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # xTB table
  output$xtb_table <- renderDT({
    df <- data.frame(
      Candidate = 1:10,
      SAScore   = c(3.125,3.160,3.200,3.239,3.241,
                    3.264,3.338,3.344,3.361,3.361),
      `ML HOMO (eV)` = c(-5.9719,-5.9719,-5.9719,
                         -5.6331,-5.6331,-5.6331,
                         -6.1831,-5.6331,-6.1831,-5.6331),
      `xTB HOMO (eV)` = c(-5.8796,-5.9140,-5.8949,
                          -6.1266,-6.0983,-6.1346,
                          -5.9269,-5.6874,-6.1366,-5.7435),
      `Error (eV)` = c(0.092,0.058,0.077,0.494,0.465,
                       0.502,0.256,0.054,0.047,0.110),
      Result = c("Pass","Pass","Pass","Fail","Fail",
                 "Fail","Fail","Pass","Pass","Pass"),
      check.names = FALSE
    )
    datatable(df,
              options = list(pageLength = 10, dom = "t"),
              rownames = FALSE) %>%
      DT::formatStyle("Result",
                      color = DT::styleEqual(c("Pass","Fail"),
                                             c("#155724","#721c24")),
                      fontWeight = "bold"
      )
  })
  
  # Hammett table
  output$hammett_table <- renderTable({
    HAMMETT[order(-HAMMETT$sigma), ]
  }, striped = TRUE, bordered = TRUE)
  
  # ── Design Engine ──────────────────────────────────────
  mol_data <- reactive({
    if (is.null(cates)) return(NULL)
    idx <- as.integer(input$mol_id)
    if (idx < 1 || idx > nrow(cates)) return(NULL)
    cates[idx]
  })
  
  output$design_output <- renderUI({
    req(input$run_btn)
    m <- mol_data()
    if (is.null(m)) {
      return(div(class="info-box",
                 "Could not load molecule data. Check BASE path in app.R"))
    }
    
    target    <- input$target_shift
    cate      <- m$CATE
    cate_se   <- m$CATE_se
    homo      <- m$homo_ev
    ewg       <- m$ewg_weighted
    mw        <- m$mol_weight
    
    req_ewg_delta <- if (cate != 0) target / cate else NA
    cyano_needed  <- if (!is.na(req_ewg_delta))
      abs(req_ewg_delta) / 0.66 else NA
    pred_homo     <- homo + cate * 0.66
    
    direction_text <- if (cate > 0)
      "POSITIVE CATE: Adding EWGs is predicted to RAISE HOMO for this scaffold."
    else
      "NEGATIVE CATE: Adding EWGs is predicted to LOWER HOMO for this scaffold."
    
    direction_class <- if (cate > 0) "warn" else "pass"
    
    tagList(
      fluidRow(
        column(3, div(class="result-box",
                      h5("Current HOMO"),
                      div(class="big-num mono", round(homo, 3)),
                      div(class="unit", "eV - experimental")
        )),
        column(3, div(class="result-box",
                      h5("Target HOMO"),
                      div(class="big-num mono",
                          style="color:#2E75B6;",
                          round(homo + target, 3)),
                      div(class="unit", "eV - desired")
        )),
        column(3, div(class="result-box",
                      h5("Per-molecule CATE"),
                      div(class="big-num mono",
                          if(cate > 0) paste0("+", round(cate, 4))
                          else round(cate, 4)),
                      div(class="unit", "eV per unit EWG score")
        )),
        column(3, div(class="result-box",
                      h5("Molecular weight"),
                      div(class="big-num mono", round(mw, 1)),
                      div(class="unit", "g/mol")
        ))
      ),
      
      div(class=paste("info-box", direction_class),
          direction_text),
      
      div(class="section-head", style="margin-top:16px;",
          "Causal Reasoning Chain"),
      
      div(class="step-box",
          div(class="step-num", "Step 1 - Per-molecule CATE"),
          p(paste0("Causal forest (4,000 trees) estimates CATE = ",
                   if(cate>0) "+" else "",
                   round(cate, 4),
                   " eV per unit EWG score for this molecule. ",
                   "This is not the population average (-0.022 eV) ",
                   "but a molecule-specific estimate. ",
                   "CATE SE = ", round(cate_se, 4), " eV."),
            style="font-size:13px; margin:0;")
      ),
      div(class="step-box",
          div(class="step-num", "Step 2 - Required EWG increase"),
          p(paste0("Required change in EWG score = target / CATE = ",
                   target, " / ", round(cate, 4), " = ",
                   if(!is.na(req_ewg_delta))
                     round(abs(req_ewg_delta), 2)
                   else "undefined (CATE = 0)"),
            style="font-size:13px; margin:0;")
      ),
      div(class="step-box",
          div(class="step-num", "Step 3 - Cyano group translation"),
          p(paste0("Each cyano group (C≡N) has Hammett sigma-para = 0.66. ",
                   "Groups needed = ",
                   if(!is.na(cyano_needed))
                     paste0(round(cyano_needed, 1),
                            " → round up to ",
                            ceiling(cyano_needed))
                   else "undefined"),
            style="font-size:13px; margin:0;")
      ),
      div(class="step-box",
          div(class="step-num", "Step 4 - Synthesizability"),
          p("Generated candidates filtered by SAScore < 4.0. All 10
           candidates from this analysis pass. SYBA score was not
           computed as the package is no longer publicly available.",
            style="font-size:13px; margin:0;")
      ),
      div(class="step-box",
          div(class="step-num", "Step 5 - xTB validation"),
          p("GFN2-xTB validates predictions. Calibration offset = 4.6197 eV
           from 8 experimental molecules. Fixed conformer seed = 42.
           MAE = 0.2155 eV across 10 candidates - 6/10 pass 0.25 eV threshold.",
            style="font-size:13px; margin:0;")
      ),
      
      div(class="section-head", style="margin-top:16px;",
          "xTB Validation - Pre-computed Top 10 Candidates"),
      DTOutput("val_table_inner")
    )
  })
  
  output$val_table_inner <- renderDT({
    df <- data.frame(
      `#` = 1:10,
      SAScore = c(3.125,3.160,3.200,3.239,3.241,
                  3.264,3.338,3.344,3.361,3.361),
      `ML HOMO` = c(-5.9719,-5.9719,-5.9719,
                    -5.6331,-5.6331,-5.6331,
                    -6.1831,-5.6331,-6.1831,-5.6331),
      `xTB HOMO` = c(-5.8796,-5.9140,-5.8949,
                     -6.1266,-6.0983,-6.1346,
                     -5.9269,-5.6874,-6.1366,-5.7435),
      `Error (eV)` = c(0.092,0.058,0.077,0.494,0.465,
                       0.502,0.256,0.054,0.047,0.110),
      Result = c("Pass","Pass","Pass","Fail","Fail",
                 "Fail","Fail","Pass","Pass","Pass"),
      check.names = FALSE
    )
    datatable(df,
              options = list(pageLength = 10, dom = "t"),
              rownames = FALSE) %>%
      DT::formatStyle("Result",
                      color = DT::styleEqual(c("Pass","Fail"),
                                             c("#155724","#721c24")),
                      fontWeight = "bold"
      )
  })
}

shinyApp(ui = ui, server = server)