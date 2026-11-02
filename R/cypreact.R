# ---------------------------------------------------------------------------
# Robust R wrapper around the CypReact command-line tool (cypreact.jar).
#
# CypReact predicts, for a given small molecule (SMILES), whether it is a
# *reactant* ("R") or a *non-reactant* ("N") of each of the nine most important
# human cytochrome-P450 (CYP450) isoforms:
#
#     1A2, 2A6, 2B6, 2C8, 2C9, 2C19, 2D6, 2E1, 3A4
#
# The underlying tool is a Java program, so this module shells out to `java`.
#
# Why a wrapper is needed
# -----------------------
# The bundled cypreact.jar has two behaviours that make naive batch use
# fragile. This wrapper works around both so result rows always line up 1:1
# with the input SMILES, in the original order:
#
#   1. A single un-parseable SMILES raises an exception that aborts the *whole*
#      batch. We detect the offending SMILES from the error message, drop it
#      (status "parse_error") and re-run until the batch succeeds.
#
#   2. Molecules that parse but are out of scope (mixtures, very large
#      molecules, sugars, lipids, ...) are *silently dropped* from the output
#      and the remaining rows are re-numbered. The tool prints
#      "Molecule with index: N is not valid ..." for each; we parse those and
#      re-insert them as status "filtered" (non-reactant for every CYP).
# ---------------------------------------------------------------------------

CYPREACT_ALL_CYPS <- c("1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4")

# Column order of the CYP predictions in the CSV produced by cypreact.jar.
.CYPREACT_OUTPUT_ORDER <- c("1A2", "2B6", "2A6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4")

# Locate the bundle (folder containing cypreact.jar and supportfiles/).
cypreact_bundle_dir <- function() {
  env <- Sys.getenv("CYPREACT_BUNDLE", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  # Default: <app dir>/cypreact_bundle. getwd() is the app directory when run
  # by shiny::runApp()/shiny-server.
  normalizePath(file.path(getwd(), "cypreact_bundle"), mustWork = FALSE)
}

cypreact_jar <- function() file.path(cypreact_bundle_dir(), "cypreact.jar")

cypreact_java <- function() {
  j <- Sys.getenv("JAVA", "")
  if (nzchar(j)) j else "java"
}

# CypReact (via CDK) uses the JNI-InChI native library, which jnati unpacks at
# runtime into "$HOME/.jnati/repo/...". On hardened hosts such as Posit Connect
# the process HOME (e.g. /home/rstudio-connect) is NOT writable, which makes the
# tool crash with "Error creating directory: .../.jnati/repo/...". To avoid any
# server configuration we point the JVM's user.home at a writable directory we
# control (reused for the whole R session so the native code is unpacked once).
cypreact_home <- function() {
  env <- Sys.getenv("CYPREACT_HOME", "")
  dir <- if (nzchar(env)) env else file.path(tempdir(), "cypreact_jnati_home")
  if (!dir.exists(dir)) {
    dir.create(file.path(dir, ".jnati", "repo"), recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(dir, mustWork = FALSE)
}

# Quick environment check; returns NULL if OK, otherwise an explanatory string.
cypreact_check_env <- function() {
  jar <- cypreact_jar()
  if (!file.exists(jar)) {
    return(paste0(
      "cypreact.jar not found at '", jar, "'. Set the CYPREACT_BUNDLE ",
      "environment variable to the folder that contains cypreact.jar and ",
      "supportfiles/."
    ))
  }
  java <- cypreact_java()
  ok <- tryCatch({
    suppressWarnings(system2(java, "-version", stdout = TRUE, stderr = TRUE))
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    return(paste0(
      "Java runtime not found (tried '", java, "'). Install a JRE (Java 8+) ",
      "or set the JAVA environment variable to the java executable."
    ))
  }
  NULL
}

# Split a block of text into a character vector of SMILES.
# Accepts one SMILES per line and/or comma-separated. Ignores blanks and
# comment (#) lines; tolerates a leading "SMILES=" and "smiles<space>name".
cypreact_parse_input <- function(text) {
  if (is.null(text) || length(text) == 0) return(character(0))
  text <- paste(text, collapse = "\n")
  text <- gsub("\r", "\n", text, fixed = TRUE)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  out <- character(0)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#")) next
    for (chunk in strsplit(line, ",", fixed = TRUE)[[1]]) {
      tok <- trimws(chunk)
      if (!nzchar(tok)) next
      if (grepl("^smiles=", tok, ignore.case = TRUE)) {
        tok <- trimws(sub("^smiles=", "", tok, ignore.case = TRUE))
      }
      # Keep only the first whitespace-separated token ("smiles name" -> smiles).
      tok <- strsplit(tok, "[[:space:]]+")[[1]][1]
      if (!is.na(tok) && nzchar(tok)) out <- c(out, tok)
    }
  }
  out
}

# Run the jar once over `smiles` for the given `cyps`. Returns a list with
# rc, output (stdout+stderr combined) and the path to the output CSV.
.cypreact_run_jar <- function(smiles, cyps, workdir) {
  in_csv  <- file.path(workdir, "input.csv")
  out_csv <- file.path(workdir, "output.csv")
  if (file.exists(out_csv)) unlink(out_csv)
  writeLines(smiles, in_csv, useBytes = TRUE)

  # JVM options: redirect user.home (and java.io.tmpdir) to writable dirs so
  # the JNI-InChI native library can be unpacked even when $HOME is read-only.
  home <- cypreact_home()
  jvm_opts <- c(
    paste0("-Duser.home=", shQuote(home)),
    paste0("-Djava.io.tmpdir=", shQuote(home))
  )
  args <- c(
    jvm_opts,
    "-jar", shQuote(cypreact_jar()),
    shQuote(paste0(cypreact_bundle_dir(), .Platform$file.sep)),
    shQuote(in_csv),
    shQuote(out_csv),
    shQuote(paste(cyps, collapse = ","))
  )
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(cypreact_bundle_dir())
  res <- suppressWarnings(
    system2(cypreact_java(), args, stdout = TRUE, stderr = TRUE)
  )
  rc <- attr(res, "status")
  if (is.null(rc)) rc <- 0L
  list(rc = rc, output = paste(res, collapse = "\n"), out_csv = out_csv)
}

# Read the output CSV into a list of rows: list(inchikey=, preds=named chr).
.cypreact_read_output <- function(out_csv) {
  if (!file.exists(out_csv)) return(list())
  df <- tryCatch(
    utils::read.csv(out_csv, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, colClasses = "character"),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(list())
  # Columns by position: 1=Inchiky, 2=SMILES, 3=Title, 4..12 = CYPs in
  # .CYPREACT_OUTPUT_ORDER. (Header names contain spaces / "2|E1", so use index.)
  rows <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    preds <- stats::setNames(
      rep(NA_character_, length(.CYPREACT_OUTPUT_ORDER)),
      .CYPREACT_OUTPUT_ORDER
    )
    for (k in seq_along(.CYPREACT_OUTPUT_ORDER)) {
      col <- 3 + k
      val <- if (col <= ncol(df)) trimws(df[i, col]) else ""
      if (nzchar(val) && tolower(val) != "null") preds[k] <- val
    }
    rows[[i]] <- list(
      inchikey = trimws(df[i, 1]),
      preds = preds
    )
  }
  rows
}

# Main entry point.
#
# Returns a data.frame with columns:
#   smiles, status, message, inchikey, <one column per requested CYP>
# status is one of "ok", "filtered", "parse_error".
cypreact_predict <- function(smiles, cyps = CYPREACT_ALL_CYPS,
                             max_parse_retries = 500L, progress = NULL) {
  smiles <- as.character(smiles)
  requested <- intersect(cyps, CYPREACT_ALL_CYPS)
  if (length(requested) == 0) requested <- CYPREACT_ALL_CYPS

  n <- length(smiles)
  make_empty <- function() {
    df <- data.frame(
      smiles = smiles,
      status = rep(NA_character_, n),
      message = rep("", n),
      inchikey = rep(NA_character_, n),
      stringsAsFactors = FALSE
    )
    for (c in requested) df[[c]] <- rep(NA_character_, n)
    df
  }
  if (n == 0) return(make_empty()[0, , drop = FALSE])

  err <- cypreact_check_env()
  if (!is.null(err)) stop(err, call. = FALSE)

  res <- vector("list", n)  # per original index -> list(status, message, inchikey, preds)

  workdir <- file.path(tempdir(), paste0("cypreact_", as.integer(Sys.time()), "_", Sys.getpid()))
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  active <- seq_len(n)              # original indices still submitted to the jar
  retries <- 0L
  repeat {
    if (length(active) == 0) break
    if (!is.null(progress)) {
      progress(sprintf("Running CypReact on %d molecule(s)...", length(active)))
    }
    submitted <- smiles[active]
    run <- .cypreact_run_jar(submitted, requested, workdir)
    combined <- run$output

    # --- Detect an un-parseable SMILES that aborted the run ---------------
    m <- regmatches(combined, regexpr("could not parse '.*?'", combined))
    aborted <- (!file.exists(run$out_csv)) || (isTRUE(run$rc != 0) && length(m) > 0)
    if (length(m) == 1 && aborted) {
      bad <- sub("^could not parse '", "", m)
      bad <- sub("'$", "", bad)
      bad <- trimws(bad)
      removed <- FALSE
      for (pos in seq_along(active)) {
        oi <- active[pos]
        if (trimws(smiles[oi]) == bad) {
          res[[oi]] <- list(
            status = "parse_error",
            message = "SMILES could not be parsed by CypReact (CDK).",
            inchikey = NA_character_, preds = NULL
          )
          active <- active[-pos]
          removed <- TRUE
          break
        }
      }
      retries <- retries + 1L
      if (removed && retries <= max_parse_retries) next
      if (!removed) {
        stop(paste0(
          "CypReact reported a SMILES parse error that could not be ",
          "attributed to an input:\n", substr(combined, 1, 2000)
        ), call. = FALSE)
      }
    }

    if (!file.exists(run$out_csv)) {
      stop(paste0(
        "CypReact failed to run (exit code ", run$rc, ").\n",
        substr(combined, 1, 2000)
      ), call. = FALSE)
    }

    # --- Success: find silently filtered molecules ------------------------
    idx_matches <- regmatches(
      combined, gregexpr("Molecule with index:\\s*([0-9]+)\\s+is not valid", combined)
    )[[1]]
    filtered_positions <- integer(0)
    if (length(idx_matches) > 0) {
      nums <- as.integer(gsub("\\D", "", idx_matches))
      filtered_positions <- nums  # 1-based positions into `submitted`
    }

    out_rows <- .cypreact_read_output(run$out_csv)
    out_ptr <- 1L
    for (pos in seq_along(active)) {
      oi <- active[pos]
      if (pos %in% filtered_positions) {
        res[[oi]] <- list(
          status = "filtered",
          message = paste0(
            "Out of CypReact's scope (e.g. mixture, very large molecule, ",
            "sugar/lipid); treated as non-reactant."
          ),
          inchikey = NA_character_,
          preds = stats::setNames(rep("N", length(requested)), requested)
        )
        next
      }
      if (out_ptr > length(out_rows)) {
        res[[oi]] <- list(
          status = "parse_error",
          message = "No output row produced by CypReact.",
          inchikey = NA_character_, preds = NULL
        )
        next
      }
      row <- out_rows[[out_ptr]]
      out_ptr <- out_ptr + 1L
      preds <- stats::setNames(
        vapply(requested, function(c) {
          v <- row$preds[[c]]
          if (is.null(v) || is.na(v)) NA_character_ else v
        }, character(1)),
        requested
      )
      res[[oi]] <- list(
        status = "ok", message = "",
        inchikey = if (nzchar(row$inchikey)) row$inchikey else NA_character_,
        preds = preds
      )
    }
    active <- integer(0)
    break
  }

  # Anything still active after the retry budget.
  for (oi in active) {
    if (is.null(res[[oi]])) {
      res[[oi]] <- list(
        status = "parse_error",
        message = "Exceeded parse-error retry budget.",
        inchikey = NA_character_, preds = NULL
      )
    }
  }

  # Assemble the data.frame.
  df <- make_empty()
  for (i in seq_len(n)) {
    r <- res[[i]]
    if (is.null(r)) next
    df$status[i]   <- r$status
    df$message[i]  <- r$message
    df$inchikey[i] <- r$inchikey
    if (!is.null(r$preds)) {
      for (c in requested) {
        if (c %in% names(r$preds)) df[[c]][i] <- r$preds[[c]]
      }
    }
  }
  df
}
