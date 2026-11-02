# =============================================================================
# CypReact Shiny app
#
# A thin, friendly web front-end around the CypReact tool. The user only needs
# to enter SMILES -- a single molecule or a whole batch -- and gets back, for
# each of the nine major human CYP450 isoforms, whether the molecule is
# predicted to be a reactant (metabolised) or not.
#
# CypReact reference:
#   Tian et al., "CypReact: A Software Tool for in Silico Reactant Prediction
#   for Human Cytochrome P450 Enzymes", J. Chem. Inf. Model. 2018.
#   Source: https://bitbucket.org/Leon_Ti/cypreact
#
# Requirements on the host: R (shiny, DT) and a Java runtime (JRE 8+).
# =============================================================================

library(shiny)
library(DT)

# Robust wrapper around cypreact.jar (defines cypreact_predict(), etc.).
source(file.path("R", "cypreact.R"), local = FALSE)

# ---- Example molecules (name -> SMILES) -------------------------------------
EXAMPLES <- c(
  "Caffeine"       = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  "Aspirin"        = "CC(=O)Oc1ccccc1C(=O)O",
  "Ibuprofen"      = "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
  "Paracetamol"    = "CC(=O)Nc1ccc(O)cc1",
  "Warfarin"       = "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
  "Diclofenac"     = "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"
)
EXAMPLE_BATCH <- paste(unname(EXAMPLES), collapse = "\n")

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background:#f6f7f9; }
    .cr-title { font-weight:700; }
    .cr-sub   { color:#5b6370; margin-top:-6px; }
    .cr-card  { background:#fff; border:1px solid #e3e6ea; border-radius:10px;
                padding:16px 18px; margin-bottom:16px; }
    .cr-R { color:#0a7d33; font-weight:700; }
    .cr-N { color:#9aa0a6; }
    .cr-badge { display:inline-block; padding:2px 8px; border-radius:10px;
                font-size:12px; font-weight:600; }
    .cr-badge-ok       { background:#e7f5ec; color:#0a7d33; }
    .cr-badge-filtered { background:#fff4e5; color:#9a6400; }
    .cr-badge-error    { background:#fdeaea; color:#b3261e; }
    .cr-legend { font-size:13px; color:#5b6370; }
    .cr-struct { display:block; margin:0 auto; }
    .cr-preview-wrap { margin:8px 0 4px 0; }
    .cr-preview { width:100%; max-width:300px; height:190px; border:1px solid #e3e6ea;
                  border-radius:8px; background:#fff; display:block; }
  ")),
  tags$script(src = "smiles-drawer.min.js"),
  tags$script(HTML("
    Shiny.addCustomMessageHandler('cypreact_draw', function(msg){
      try {
        var el = document.getElementById(msg.id);
        if (!el || typeof SmiDrawer === 'undefined') return;
        if (!msg.smiles) { el.innerHTML = ''; return; }
        var d = new SmiDrawer({width: msg.width || 300, height: msg.height || 190}, {});
        d.draw(msg.smiles, '#' + msg.id, msg.theme || 'light', null, function(err){
          el.innerHTML = '<text x=\"8\" y=\"20\" fill=\"#b3261e\" font-size=\"12\">invalid SMILES</text>';
        });
      } catch(e) {}
    });
    // Re-draw all [data-smiles] cells after a DataTables redraw (paging/sort).
    window.cypreactDrawTable = function(){
      try { if (typeof SmilesDrawer !== 'undefined') {
        SmilesDrawer.apply({width:120, height:84}, {}, 'data-smiles', 'light');
      } } catch(e) {}
    };
  "))),

  titlePanel(
    title = div(
      div(class = "cr-title", "CypReact – CYP450 Reactant Prediction"),
      div(class = "cr-sub",
          "Enter SMILES to predict which human cytochrome-P450 enzymes ",
          "metabolise your molecule(s).")
    ),
    windowTitle = "CypReact"
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      radioButtons(
        "mode", "Input mode",
        choices = c("Single molecule" = "single", "Batch" = "batch"),
        selected = "single", inline = TRUE
      ),

      conditionalPanel(
        "input.mode == 'single'",
        textInput("smiles_single", "SMILES",
                  value = unname(EXAMPLES[["Caffeine"]]),
                  placeholder = "e.g. CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        div(style = "margin-bottom:10px;",
            tags$small("Examples: "),
            lapply(names(EXAMPLES), function(nm) {
              actionLink(paste0("ex_", make.names(nm)), nm,
                         style = "margin-right:8px;")
            })
        ),
        tags$small("Structure preview"),
        div(class = "cr-preview-wrap",
            HTML('<svg id="cr-preview" class="cr-preview"></svg>'))
      ),

      conditionalPanel(
        "input.mode == 'batch'",
        textAreaInput("smiles_batch", "SMILES (one per line, or comma-separated)",
                      value = EXAMPLE_BATCH, rows = 8,
                      placeholder = "CCO\nc1ccccc1\n..."),
        fileInput("smiles_file", "...or upload a .csv / .txt / .smi file",
                  accept = c(".csv", ".txt", ".smi")),
        tags$small("Uploaded files are appended to the text above. One SMILES per line.")
      ),

      hr(),
      checkboxGroupInput(
        "cyps", "CYP450 isoforms to predict",
        choices = CYPREACT_ALL_CYPS, selected = CYPREACT_ALL_CYPS, inline = TRUE
      ),
      div(style = "margin:-6px 0 10px 0;",
          actionLink("cyp_all", "all"), " / ",
          actionLink("cyp_none", "none")),

      actionButton("run", "Predict", class = "btn-primary", width = "100%",
                   icon = icon("flask")),
      br(), br(),
      div(class = "cr-legend",
          strong("R"), " = predicted reactant (metabolised)  ·  ",
          strong("N"), " = non-reactant")
    ),

    mainPanel(
      width = 8,
      uiOutput("summary"),
      div(class = "cr-card",
          div(style = "display:flex; justify-content:space-between; align-items:center;",
              h4("Results", style = "margin:0;"),
              downloadButton("download", "Download CSV", class = "btn-sm")),
          br(),
          DT::dataTableOutput("results"),
          br(),
          div(class = "cr-legend",
              strong("Legend"), br(),
              span(class = "cr-R", "R"),
              " = predicted reactant (metabolised by that CYP450 isoform)",
              span(style = "margin:0 10px;", "·"),
              span(class = "cr-N", "N"),
              " = non-reactant", br(),
              span(class = "cr-badge cr-badge-ok", "predicted"),
              " molecule processed normally",
              span(style = "margin:0 10px;", "·"),
              span(class = "cr-badge cr-badge-filtered", "out of scope"),
              " e.g. mixture / very large / sugar / lipid – treated as non-reactant",
              span(style = "margin:0 10px;", "·"),
              span(class = "cr-badge cr-badge-error", "parse error"),
              " SMILES could not be read",
              span(style = "margin:0 10px;", "·"),
              span(class = "cr-badge cr-badge-error", "error"),
              " could not be predicted – skipped, batch continues"
          )
      ),
      div(class = "cr-legend",
          "Powered by CypReact (Tian et al., 2018). Predictions are ",
          "in-silico estimates and not a substitute for experimental data.")
    )
  )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  # Example quick-fill links (single mode).
  lapply(names(EXAMPLES), function(nm) {
    observeEvent(input[[paste0("ex_", make.names(nm))]], {
      updateTextInput(session, "smiles_single", value = unname(EXAMPLES[[nm]]))
    })
  })

  observeEvent(input$cyp_all,  updateCheckboxGroupInput(session, "cyps", selected = CYPREACT_ALL_CYPS))
  observeEvent(input$cyp_none, updateCheckboxGroupInput(session, "cyps", selected = character(0)))

  # Live 2D structure preview for the single-molecule SMILES (client-side).
  observe({
    input$mode  # redraw when switching back to single mode
    smi <- cypreact_parse_input(input$smiles_single)
    session$sendCustomMessage("cypreact_draw", list(
      id = "cr-preview",
      smiles = if (length(smi) > 0) smi[1] else "",
      width = 300, height = 190, theme = "light"
    ))
  })

  # Append uploaded file contents to the batch text area.
  observeEvent(input$smiles_file, {
    f <- input$smiles_file
    if (is.null(f)) return()
    txt <- tryCatch(paste(readLines(f$datapath, warn = FALSE), collapse = "\n"),
                    error = function(e) "")
    cur <- input$smiles_batch
    new <- if (nzchar(trimws(cur))) paste(cur, txt, sep = "\n") else txt
    updateTextAreaInput(session, "smiles_batch", value = new)
  })

  results <- eventReactive(input$run, {
    smiles <- if (identical(input$mode, "single")) {
      cypreact_parse_input(input$smiles_single)
    } else {
      cypreact_parse_input(input$smiles_batch)
    }

    validate(need(length(smiles) > 0, "Please enter at least one SMILES string."))
    validate(need(length(input$cyps) > 0, "Please select at least one CYP450 isoform."))

    withProgress(message = "Running CypReact...", value = 0, {
      tryCatch(
        cypreact_predict(
          smiles, cyps = input$cyps,
          progress = function(done, total, label) {
            setProgress(value = if (total > 0) done / total else 0, detail = label)
          }
        ),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = NULL)
          NULL
        }
      )
    })
  })

  output$summary <- renderUI({
    df <- results()
    if (is.null(df)) return(NULL)
    n       <- nrow(df)
    n_ok    <- sum(df$status == "ok")
    n_filt  <- sum(df$status == "filtered")
    n_err   <- sum(df$status %in% c("parse_error", "error"))
    div(class = "cr-card",
        span(class = "cr-badge cr-badge-ok",
             sprintf("%d predicted", n_ok)),
        if (n_filt > 0) span(class = "cr-badge cr-badge-filtered",
             style = "margin-left:8px;", sprintf("%d out of scope", n_filt)),
        if (n_err > 0) span(class = "cr-badge cr-badge-error",
             style = "margin-left:8px;", sprintf("%d could not be processed", n_err)),
        span(style = "margin-left:10px; color:#5b6370;",
             sprintf("· %d molecule(s) total", n))
    )
  })

  # Build the display table.
  display_df <- reactive({
    df <- results()
    if (is.null(df)) return(NULL)
    cyp_cols <- intersect(CYPREACT_ALL_CYPS, names(df))
    status_label <- c(ok = "predicted", filtered = "out of scope",
                      parse_error = "parse error", error = "error")[df$status]
    out <- data.frame(
      "#" = seq_len(nrow(df)),
      SMILES = df$smiles,
      Status = unname(status_label),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    for (c in cyp_cols) out[[c]] <- df[[c]]
    out$InChIKey <- df$inchikey
    out$Note <- df$message
    out
  })

  # Escape a string so it is safe inside an HTML attribute.
  esc_attr <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    x <- gsub(">", "&gt;",  x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
  }

  output$results <- DT::renderDataTable({
    out <- display_df()
    if (is.null(out)) return(NULL)
    cyp_cols <- intersect(CYPREACT_ALL_CYPS, names(out))

    # Insert a client-side 2D structure column (drawn by smiles-drawer).
    struct <- vapply(out$SMILES, function(s) {
      sprintf('<svg class="cr-struct" data-smiles="%s" width="120" height="84"></svg>',
              esc_attr(s))
    }, character(1), USE.NAMES = FALSE)
    out$Structure <- struct
    ordered <- c("#", "SMILES", "Structure", "Status", cyp_cols, "InChIKey", "Note")
    out <- out[, intersect(ordered, names(out)), drop = FALSE]

    # We render raw HTML (the SVG column), so escape is FALSE globally; the only
    # user-supplied text column (SMILES) is HTML-escaped by hand here.
    out$SMILES <- esc_attr(out$SMILES)

    # 0-based column indices to centre (structure + CYP columns).
    center0 <- match(c("Structure", cyp_cols), names(out)) - 1L

    dt <- DT::datatable(
      out, rownames = FALSE, escape = FALSE,
      options = list(
        pageLength = 25, scrollX = TRUE,
        columnDefs = list(list(className = "dt-center", targets = center0)),
        # Re-render all structure SVGs after each DataTables draw (paging/sort).
        drawCallback = DT::JS("function(settings){ if(window.cypreactDrawTable) window.cypreactDrawTable(); }")
      )
    )
    # Colour R / N cells.
    for (c in cyp_cols) {
      dt <- DT::formatStyle(
        dt, c,
        color = DT::styleEqual(c("R", "N"), c("#0a7d33", "#9aa0a6")),
        fontWeight = DT::styleEqual(c("R", "N"), c("bold", "normal")),
        backgroundColor = DT::styleEqual(c("R"), c("#e7f5ec"))
      )
    }
    # Tint the status of failed / out-of-scope rows.
    dt <- DT::formatStyle(
      dt, "Status",
      color = DT::styleEqual(
        c("predicted", "out of scope", "parse error", "error"),
        c("#0a7d33", "#9a6400", "#b3261e", "#b3261e")
      ),
      fontWeight = "bold"
    )
    dt
  })

  output$download <- downloadHandler(
    filename = function() paste0("cypreact_results_", Sys.Date(), ".csv"),
    content = function(file) {
      out <- display_df()
      if (is.null(out)) out <- data.frame()
      utils::write.csv(out, file, row.names = FALSE, na = "")
    }
  )
}

shinyApp(ui, server)
