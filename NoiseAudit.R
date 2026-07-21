##################################################
### AVIATION FORECAST NOISE ANALYSIS
########################################################
########################################################
########################################################
library(data.table)
library(stringr)
library(lubridate)
library(ggplot2)

#########
# FUNCTION 1 — Parse TAF (INCLUDING TEMPO/INTER)
expand_taf <- function(file){

taf <- fread(file, header=FALSE, sep="\n")
out <- data.table()

for(i in 1:nrow(taf)){
  
  line <- taf$V1[i]
  
  date_str <- str_extract(line, "^\\d{8}")
  time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
  
  issued <- ymd_hms(paste(date_str, time_str))
  if(is.na(issued)) next
  
  # ✅ FIXED HAZARD DETECTION
  has_fg <- str_detect(line,"\\bFG\\b")
  has_ts <- str_detect(line,"\\bTS\\b|TSRA|VCTS|TSGR")
  
  vis_groups <- str_extract_all(
    line,
    "(?<!\\d)\\d{4}(?!\\d)"
  )[[1]]
  
  vis_groups <- suppressWarnings(as.numeric(vis_groups))
  
  has_vis <- any(
    vis_groups <= 5000 &
      vis_groups >= 0,
    na.rm = TRUE
  )
  
  cld_vals <- str_extract_all(
    line,
    "\\b(?:BKN|OVC)\\d{3}\\b"
  )[[1]]
  
  cld_vals <- as.numeric(
    str_extract(cld_vals, "\\d{3}")
  ) * 100
  
  has_cld1000 <- any(cld_vals < 1000, na.rm=TRUE)
  
  times <- seq(
    floor_date(issued, "hour"),
    floor_date(issued + hours(12), "hour"),
    by="1 hour"
  )
  
  dt <- data.table(
    Issued = issued,
    validTime = times,
    fg_fcst = has_fg,
    ts_fcst = has_ts,
    vis_fcst = has_vis,
    cld1000_fcst = has_cld1000
  )
  
  out <- rbind(out, dt)
}

return(out)
}

# FUNCTION 2 — Parse METAR → events
parse_metar <- function(file){
  
  met <- fread(file, header=FALSE, sep="\n")
  
  dt <- met[, {
    
    line <- V1
    
    # Extract datetime
    date_str <- str_extract(line, "^\\d{8}")
    time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    datetime <- suppressWarnings(ymd_hms(paste(date_str, time_str)))
    
    # If bad parse, return NA row (NOT return())
    tokens <- str_split(line, " ")[[1]]
    tokens <- tokens[tokens != "" & tokens != "////"]
    
    vis <- str_extract(line, "(?<=VIS:)\\d{4}")
    vis <- suppressWarnings(as.numeric(vis))
    
    wx <- paste(tokens[grep("FG|TS|RA|SH|DZ", tokens)], collapse=" ")
    
    cld <- tokens[grep("FEW|SCT|BKN|OVC", tokens)]
    cld_base <- suppressWarnings(
      min(as.numeric(str_sub(cld,4,6))*100, na.rm=TRUE)
    )
    
    # always return SAME structure
    list(
      time = datetime,
      vis = vis,
      wx = wx,
      cld_base = cld_base
    )
    
  }, by=1:nrow(met)]
  
  # Clean bad rows
  dt <- dt[!is.na(time)]
  
  # Align to hourly grid
  dt[, time := floor_date(time, "hour")]
  
  # Hazard definitions (Centralised)
  dt[, `:=`(
    event_fg = str_detect(wx,"FG"),
    event_ts = str_detect(wx,"TS"),
    event_vis_lt1000 = !is.na(vis) & vis < 1000,
    event_cld_lt1000 = !is.infinite(cld_base) & cld_base < 1000,
    event_cld_lt2500 = !is.infinite(cld_base) & cld_base < 2500
  )]
  
  # Define event
  dt[, event_any := event_fg | event_ts | event_vis_lt1000 | event_cld_lt1000]
  
  return(dt[event_any == TRUE])
  
}

# FUNCTION 3 — Compute lead time
compute_lead_times <- function(taf_dt, met_dt){
  
  results <- data.table()
  
  # collapse TAF to issue-level
  taf_summary <- taf_dt[, .(
    fg = any(fg_fcst),
    ts = any(ts_fcst),
    vis = any(vis_fcst),
    cld1000 = any(cld1000_fcst)
  ), by = Issued]
  
  for(i in 1:nrow(met_dt)){
    
    t_event <- met_dt$time[i]
    
    # identify hazard type for this METAR
    hazard <- if(met_dt$event_fg[i]) "FG" else
      if(met_dt$event_ts[i]) "TS" else
        if(met_dt$event_vis_lt1000[i]) "VIS" else
          if(met_dt$event_cld_lt1000[i]) "CLD" else NA
    
    if(is.na(hazard)) next
    
    # select matching forecasts
    candidates <- taf_summary[
      Issued <= t_event &
        Issued >= t_event - hours(24) &
        (
          (hazard=="FG"  & fg) |
            (hazard=="TS"  & ts) |
            (hazard=="VIS" & vis) |
            (hazard=="CLD" & cld1000)
        )
    ]
    
    if(nrow(candidates)==0) next
    
    first <- candidates[which.min(Issued)]
    last  <- candidates[which.max(Issued)]
    
    results <- rbind(results, data.table(
      event_time = t_event,
      hazard = hazard,
      first_issue = first$Issued,
      last_issue  = last$Issued,
      lead_first = as.numeric(difftime(t_event, first$Issued, units="hours")),
      lead_last  = as.numeric(difftime(t_event, last$Issued, units="hours")),
      signal_gap = as.numeric(difftime(last$Issued, first$Issued, units="hours"))
    ))
  }
  
  return(results)
}

# FUNCTION 4 - Run Full Analysis
run_analysis <- function(taf_files, metar_files){
  
  all_results <- data.table()
  
  for(k in seq_along(taf_files)){
    
    cat("Processing:", taf_files[k], "\n")
    
    taf_dt <- expand_taf(taf_files[k])
    met_dt <- parse_metar(metar_files[k])
    
    res <- compute_lead_times(taf_dt, met_dt)
    
    if(nrow(res)==0) next
    
    res[, aerodrome := toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    #res[, aerodrome := toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res)
  }
  
  return(all_results)
}

# FUNCTION 5 - EXPAND TAF WITH AMD
expand_taf_with_amd <- function(file){

taf <- fread(file, header=FALSE, sep="\n")
out <- data.table()

for(i in 1:nrow(taf)){
  
  line <- taf$V1[i]
  
  date_str <- str_extract(line, "^\\d{8}")
  time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
  
  issued <- ymd_hms(paste(date_str, time_str))
  if(is.na(issued)) next
  
  is_amd <- str_detect(line, "\\bAMD\\b")
  
  # ✅ FIXED TS
  has_fg <- str_detect(line,"\\bFG\\b")
  has_ts <- str_detect(line,"\\bTS\\b|TSRA|VCTS|TSGR")
  
  vis_groups <- str_extract_all(
    line,
    "(?<!\\d)\\d{4}(?!\\d)"
  )[[1]]
  
  vis_groups <- suppressWarnings(as.numeric(vis_groups))
  
  has_vis <- any(
    vis_groups <= 1000 &
      vis_groups >= 0,
    na.rm = TRUE
  )
  
  cld_vals <- str_extract_all(
    line,
    "\\b(?:BKN|OVC)\\d{3}\\b"
  )[[1]]
  
  cld_vals <- as.numeric(
    str_extract(cld_vals, "\\d{3}")
  ) * 100
  
  has_cld1000 <- any(cld_vals < 1000, na.rm=TRUE)
  
  times <- seq(
    floor_date(issued, "hour"),
    floor_date(issued + hours(12), "hour"),
    by="1 hour"
  )
  
  dt <- data.table(
    Issued = issued,
    validTime = times,
    is_amd = is_amd,
    fg_fcst = has_fg,
    ts_fcst = has_ts,
    vis_fcst = has_vis,
    cld1000_fcst = has_cld1000
  )
  
  out <- rbind(out, dt, fill=TRUE)
}

return(out)
}

# FUNCTION 6 - AMD-AWARE LEAD TIME CALCULATION
compute_amd_analysis <- function(taf_dt, met_dt){
  
  results <- data.table()
  
  for(i in 1:nrow(met_dt)){
    
    t_event <- met_dt$time[i]
    
    hazard <- if(met_dt$event_fg[i]) "FG" else
      if(met_dt$event_ts[i]) "TS" else
        if(met_dt$event_vis_lt1000[i]) "VIS" else
          if(met_dt$event_cld_lt1000[i]) "CLD" else NA
    
    if(is.na(hazard)) next
    
    # CRITICAL: match forecasts at EVENT TIME
    candidates <- taf_dt[
      validTime >= t_event - hours(1) &
        validTime <= t_event + hours(1) &
        Issued <= t_event &
        Issued >= t_event - hours(24)
    ]
    
    candidates <- candidates[
      order(Issued, validTime),
      .SD[.N],
      by = Issued
    ]
    
    if(nrow(candidates) < 2) next
    
    candidates <- candidates[order(Issued)]
    
    before <- candidates[.N - 1]
    after  <- candidates[.N]
    
    # forecast presence BEFORE (at event time)
    before_fcst <- (
      (hazard=="FG"  & before$fg_fcst) |
        (hazard=="TS"  & before$ts_fcst) |
        (hazard=="VIS" & before$vis_fcst) |
        (hazard=="CLD" & before$cld1000_fcst)
    )
    
    # forecast presence AFTER
    after_fcst <- (
      (hazard=="FG"  & after$fg_fcst) |
        (hazard=="TS"  & after$ts_fcst) |
        (hazard=="VIS" & after$vis_fcst) |
        (hazard=="CLD" & after$cld1000_fcst)
    )
    
    before_any <- before$fg_fcst | before$ts_fcst | before$vis_fcst | before$cld1000_fcst
    after_any  <- after$fg_fcst  | after$ts_fcst  | after$vis_fcst  | after$cld1000_fcst
    
    # timing distance to event
    time_diff_before <- abs(difftime(t_event, before$Issued, units="hours"))
    time_diff_after  <- abs(difftime(t_event, after$Issued, units="hours"))
    
    timing_improved <- time_diff_after < time_diff_before
    
    added_signal <- (!before_any & after_any)
    removed_signal <- (before_any & !after_any)
    
    # relaxed improvement (includes timing refinement)
    useful_update <- (
      added_signal |
        timing_improved |
        (!before_fcst & after_fcst)
    )
    
    improved  <- (!before_fcst & after_fcst)
    degraded  <- (before_fcst & !after_fcst)
    unchanged <- (before_fcst == after_fcst)
    
    results <- rbind(results, data.table(
      event_time = t_event,
      hazard = hazard,
      taf_issue = after$Issued,
      lead_time = as.numeric(difftime(t_event, after$Issued, units="hours")),
      is_amd = after$is_amd,
      useful_update = useful_update,
      improved = improved,
      degraded = degraded,
      unchanged = unchanged
    ))
  }
  
  return(results)
}

# FUNCTION 7 - RUN AMD ANALYSIS
run_amd_analysis <- function(taf_files, metar_files){
  
  all_results <- data.table()
  
  for(k in seq_along(taf_files)){
    
    cat("Processing AMD:", taf_files[k], "\n")
    
    taf_dt <- expand_taf_operational(taf_files[k])
    met_dt <- parse_metar(metar_files[k])
    
    res <- compute_amd_analysis(taf_dt, met_dt)
    
    if(nrow(res)==0) next
    
    res[, aerodrome := toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    #res[, aerodrome := toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res)
  }
  
  return(all_results)
}

# FUNCTION 8 - Extract TEMPO / INTER / PROB forecasts
expand_taf_prob_tempo <- function(file){

taf <- fread(file, header = FALSE, sep = "\n")
out <- data.table()

for(i in 1:nrow(taf)){
  
  line <- taf$V1[i]
  
  date_str <- str_extract(line, "^\\d{8}")
  time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
  issued <- ymd_hms(paste(date_str, time_str))
  if(is.na(issued)) next
  
  is_tempo <- str_detect(line, "\\bTEMPO\\b") | str_detect(line, "\\bINTER\\b")
  is_prob30 <- str_detect(line, "PROB30")
  is_prob40 <- str_detect(line, "PROB40")
  
  # ✅ FIXED TS
  has_fg <- str_detect(line,"\\bFG\\b")
  has_ts <- str_detect(line,"\\bTS\\b|TSRA|VCTS|TSGR")
  
  vis_groups <- str_extract_all(
    line,
    "(?<!\\d)\\d{4}(?!\\d)"
  )[[1]]
  
  vis_groups <- suppressWarnings(as.numeric(vis_groups))
  
  has_vis <- any(
    vis_groups <= 1000 &
      vis_groups >= 0,
    na.rm = TRUE
  )
  
  cld_vals <- str_extract_all(
    line,
    "\\b(?:BKN|OVC)\\d{3}\\b"
  )[[1]]
  
  cld_vals <- as.numeric(
    str_extract(cld_vals, "\\d{3}")
  ) * 100
  has_cld <- any(cld_vals < 1000, na.rm=TRUE)
  
  times <- seq(
    floor_date(issued, "hour"),
    floor_date(issued + hours(12), "hour"),
    by = "1 hour"
  )
  
  dt <- data.table(
    Issued = issued,
    validTime = times,
    is_tempo = is_tempo,
    is_prob30 = is_prob30,
    is_prob40 = is_prob40,
    fg_fcst = has_fg,
    ts_fcst = has_ts,
    vis_fcst = has_vis,
    cld_fcst = has_cld
  )
  
  out <- rbind(out, dt, fill=TRUE)
}

return(out)
}

# FUNCTION 9 - When TEMP/INTER forecast exists → does hazard occur?
verify_tempo <- function(taf_dt, met_dt){
  
  results <- data.table()
  
  # use ALL times
  met_lookup <- met_dt[, .(
    event_any = any(event_any)
  ), by=time]
  
  for(i in 1:nrow(taf_dt)){
    
    if(!taf_dt$is_tempo[i]) next
    
    t <- taf_dt$validTime[i]
    
    obs_row <- met_lookup[time == t]
    if(nrow(obs_row)==0) next
    
    fcst_any <- taf_dt$fg_fcst[i] | taf_dt$ts_fcst[i] |
      taf_dt$vis_fcst[i] | taf_dt$cld_fcst[i]
    
    obs_any <- obs_row$event_any
    
    results <- rbind(results, data.table(
      forecast = fcst_any,
      observed = obs_any,
      hit = fcst_any & obs_any,
      false_alarm = fcst_any & !obs_any,
      miss = !fcst_any & obs_any
    ))
  }
  
  return(results)
}

# FUNCTION 10 - Verify PROB forecasts
# If PROB30 used, does event occur ~30% of time?
verify_prob <- function(taf_dt, met_dt){
  
  results <- data.table(
    prob = numeric(),
    observed = logical()
  )
  
  # Filter only rows where PROB exists AND hazard is forecast
  taf_dt <- taf_dt[
    !is.na(prob_type) &
      (fg_fcst | ts_fcst | vis_fcst | cld_fcst)
  ]
  
  met_lookup <- met_dt[, .(
    event_any = any(event_any)
  ), by=time]
  
  for(i in 1:nrow(taf_dt)){
    
    t <- taf_dt$validTime[i]
    
    obs_row <- met_lookup[time == t]
    if(nrow(obs_row) == 0) next
    
    obs_any <- obs_row$event_any
    
    results <- rbind(results, data.table(
      prob = taf_dt$prob_type[i],
      observed = obs_any
    ))
  }
  
  return(results)
}

# FUNCTION 11 - TEMPO / INTER Run across aerodromes
run_tempo_analysis <- function(taf_files, metar_files){
  
  all_results <- data.table()
  
  for(k in seq_along(taf_files)){
    
    tempo_events <- extract_tempo_events(taf_files[k])
    met_dt <- parse_metar_full(metar_files[k])
    
    res <- verify_tempo_events(tempo_events, met_dt)
    
    if(nrow(res)==0) next
    
    res[, aerodrome :=
          toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res, fill=TRUE)
  }
  
  return(all_results)
}

# FUNCTION 12 - PROB Run across aerodromes
run_prob_analysis <- function(taf_files, metar_files){
  
  all_results <- data.table()
  
  for(k in seq_along(taf_files)){
    
    prob_events <- extract_prob_events(taf_files[k])
    met_dt <- parse_metar_full(metar_files[k])
    
    res <- verify_prob_events(prob_events, met_dt)
    
    if(nrow(res) == 0) next
    
    res[, aerodrome :=
          toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res, fill=TRUE)
  }
  
  return(all_results)
}

# FUNCTION 13 - Create full METAR dataset
parse_metar_full <- function(file){
  
  met <- fread(
    file,
    header = FALSE,
    sep = "\n"
  )
  
  dt <- met[, {
    
    line <- V1
    
    date_str <- str_extract(
      line,
      "^\\d{8}"
    )
    
    time4 <- str_match(
      line,
      "^\\d{8}\\s+(\\d{4})"
    )[,2]
    
    datetime <- as.POSIXct(
      paste0(
        substr(date_str,1,4), "-",
        substr(date_str,5,6), "-",
        substr(date_str,7,8), " ",
        substr(time4,1,2), ":",
        substr(time4,3,4), ":00"
      ),
      tz = "UTC"
    )
    
    tokens <- str_split(line, " ")[[1]]
    tokens <- tokens[
      tokens != "" &
        tokens != "////"
    ]
    
    wind_match <- str_match(
      line,
      "\\b(\\d{3})(\\d{2})/(\\d{2})KT\\b"
    )
    
    wind_dir <- suppressWarnings(
      as.numeric(wind_match[2])
    )
    
    wind_spd <- suppressWarnings(
      as.numeric(wind_match[3])
    )
    
    temp <- suppressWarnings(
      as.numeric(
        str_extract(
          line,
          "-?\\d+\\.\\d+(?=/)"
        )
      )
    )
    
    vis <- str_extract(
      line,
      "(?<=VIS:)\\d{4}"
    )
    
    vis <- suppressWarnings(
      as.numeric(vis)
    )
    
    wx <- paste(
      tokens[
        grep(
          "FG|TS|RA|SH|DZ",
          tokens
        )
      ],
      collapse = " "
    )
    
    cld <- tokens[
      grep(
        "FEW|SCT|BKN|OVC",
        tokens
      )
    ]
    
    cld_base <- suppressWarnings(
      min(
        as.numeric(
          str_sub(
            cld,
            4,
            6
          )
        ) * 100,
        na.rm = TRUE
      )
    )
    
    cloud_info <- extract_cloud_layers(line)
    
    wx_info <- extract_weather_phenomena(line)
    
    list(
      time = datetime,
      vis = vis,
      wx = wx,
      cld_base = cld_base,
      
      wind_dir = wind_dir,
      wind_spd = wind_spd,
      temp = temp,
      
      cloud_amount =
        cloud_info$cloud_amount,
      
      cloud_base =
        cloud_info$cloud_base,
      
      event_ra = wx_info$ra,
      event_shra = wx_info$shra,
      event_dz = wx_info$dz,
      event_br = wx_info$br,
      event_gr = wx_info$gr,
      event_gs = wx_info$gs
    )
    
  }, by = 1:nrow(met)]
  
  dt <- dt[!is.na(time)]
  
  dt[, time := floor_date(
    time,
    "hour"
  )]
  
  dt[, `:=`(
    
    event_fg =
      str_detect(
        wx,
        "FG|BCFG|MIFG|PRFG"
      ),
    
    event_ts =
      str_detect(
        wx,
        "\\bTS\\b|TSRA|VCTS|TSGR"
      ),
    
    event_vis_lt1000 =
      !is.na(vis) &
      vis < 1000,
    
    event_cld_lt1000 =
      !is.infinite(cld_base) &
      cld_base < 1000
    
  )]
  
  dt[, event_any :=
       event_fg |
       event_ts |
       event_vis_lt1000 |
       event_cld_lt1000]
  
  dt
}

# FUNCTION 14 - Extract TEMPO time windows
extract_tempo_windows <- function(line){
  
  matches <- str_match_all(line, "TEMPO\\s+(\\d{4})/(\\d{4})")[[1]]
  
  if(nrow(matches) == 0) return(NULL)
  
  data.table(
    start = matches[,2],
    end   = matches[,3]
  )
}

# FUNCTION 15 - Extract PROB time windows
extract_prob_windows <- function(line){
  
  matches <- str_match_all(line, "PROB(30|40)\\s+(\\d{4})/(\\d{4})")[[1]]
  
  if(nrow(matches) == 0) return(NULL)
  
  data.table(
    prob  = as.numeric(matches[,2]),
    start = matches[,3],
    end   = matches[,4]
  )
}

# FUNCTION 16 - Build time-aware TEMPO expansion
# To replace earlier TEMPO logic
expand_taf_tempo_precise <- function(file){
  
  taf <- fread(file, header=FALSE, sep="\n")
  out <- data.table()
  
  for(i in 1:nrow(taf)){
    
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if(is.na(issued)) next
    
    # base hourly timeline
    times <- seq(
      floor_date(issued, "hour"),
      floor_date(issued + hours(12), "hour"),
      by="1 hour"
    )
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      hour = hour(times),   # ✅ ADD THIS LINE
      tempo_active = FALSE,
      fg_fcst = FALSE,
      ts_fcst = FALSE,
      vis_fcst = FALSE,
      cld_fcst = FALSE
    )
    
    # check TEMPO window
    tw <- extract_tempo_window(line)
    
    if(!is.null(tw)){
      
      # approximate conversion (simplified UTC day handling)
      start_hour <- as.numeric(substr(tw$start,3,4))
      end_hour   <- as.numeric(substr(tw$end,3,4))
      
      dt[hour >= start_hour & hour <= end_hour, tempo_active := TRUE]
      
      # hazards ONLY applied in TEMPO window
      dt[tempo_active == TRUE, `:=`(
        fg_fcst = str_detect(line,"\\bFG\\b"),
        ts_fcst = str_detect(line,"\\bTS\\b"),
        vis_fcst = any(as.numeric(str_extract_all(line,"\\b\\d{4}\\b")[[1]]) < 1000, na.rm=TRUE),
        cld_fcst = any(as.numeric(str_extract_all(line,"\\d{3}")[[1]])*100 < 1000, na.rm=TRUE)
      )]
    }
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  return(out)
}

# FUNCTION 17 - Build time-aware PROB expansion
expand_taf_prob_precise <- function(file){
  
  taf <- fread(file, header=FALSE, sep="\n")
  out <- data.table()
  
  for(i in 1:nrow(taf)){
    
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if(is.na(issued)) next
    
    times <- seq(
      floor_date(issued, "hour"),
      floor_date(issued + hours(12), "hour"),
      by="1 hour"
    )
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      prob_type = NA_real_,
      fg_fcst = FALSE,
      ts_fcst = FALSE,
      vis_fcst = FALSE,
      cld_fcst = FALSE
    )
    
    pw <- extract_prob_window(line)
    
    if(!is.null(pw)){
      
      start_hour <- as.numeric(substr(pw$start,3,4))
      end_hour   <- as.numeric(substr(pw$end,3,4))
      
      dt[, hour := hour(validTime)]
      
      dt[hour >= start_hour & hour <= end_hour, prob_type := pw$prob]
      
      # apply hazard ONLY in window
      dt[!is.na(prob_type), `:=`(
        fg_fcst = str_detect(line,"\\bFG\\b"),
        ts_fcst = str_detect(line,"\\bTS\\b"),
        vis_fcst = any(as.numeric(str_extract_all(line,"\\b\\d{4}\\b")[[1]]) < 1000, na.rm=TRUE),
        cld_fcst = any(as.numeric(str_extract_all(line,"\\d{3}")[[1]])*100 < 1000, na.rm=TRUE)
      )]
    }
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  return(out)
}

# FUNCTION 18 - Verify TEMPO precise
verify_tempo_precise <- function(taf_dt, met_dt){
  
  met_lookup <- met_dt[, .(
    event_any = any(event_any)
  ), by=time]
  
  results <- merge(
    taf_dt,
    met_lookup,
    by.x="validTime",
    by.y="time",
    all.x=TRUE
  )
  
  results[, event_any := ifelse(is.na(event_any), FALSE, event_any)]
  
  results[, forecast := fg_fcst | ts_fcst | vis_fcst | cld_fcst]
  
  results[, `:=`(
    hit = forecast & event_any,
    false_alarm = forecast & !event_any,
    miss = !forecast & event_any
  )]
  
  return(results)
}

# FUNCTION 19 - collapse TEMPO to window level
verify_tempo_events <- function(tempo_events, met_dt){
  
  results <- data.table(
    forecast = logical(),
    observed = logical()
  )
  
  for(i in 1:nrow(tempo_events)){
    
    start <- tempo_events$start_time[i]
    end   <- tempo_events$end_time[i]
    
    obs <- met_dt[time >= start & time <= end]
    
    obs_any <- if(nrow(obs)==0) FALSE else any(obs$event_any, na.rm=TRUE)
    
    results <- rbind(results, data.table(
      forecast = tempo_events$forecast[i],
      observed = obs_any
    ))
  }
  
  results[, `:=`(
    hit = forecast & observed,
    false_alarm = forecast & !observed,
    miss = !forecast & observed
  )]
  
  return(results)
}

# FUNCTION 20 - collapse TEMPO to window level
verify_prob_window <- function(taf_dt, met_dt){
  
  results <- data.table(
    prob = numeric(),
    observed = logical()
  )
  
  # ✅ Identify UNIQUE PROB windows
  prob_windows <- taf_dt[
    !is.na(prob_type)
  ][, .(
    start_time = min(validTime),
    end_time   = max(validTime),
    prob       = unique(prob_type)[1]
  ), by = .(Issued, prob_type)]
  
  # ✅ IMPORTANT: grouping by Issued AND prob_type
  # ensures PROB30 and PROB40 handled separately
  # but still not enough if multiple windows exist (see below)
  
  for(i in 1:nrow(prob_windows)){
    
    start <- prob_windows$start_time[i]
    end   <- prob_windows$end_time[i]
    
    obs <- met_dt[time >= start & time <= end]
    
    obs_any <- if(nrow(obs) == 0) FALSE else any(obs$event_any, na.rm=TRUE)
    
    results <- rbind(results, data.table(
      prob = prob_windows$prob[i],
      observed = obs_any
    ))
  }
  
  return(results)
}

# FUNCTION 21 - Extract individual PROB events
extract_prob_events <- function(file){
  
  taf_blocks <- build_taf_blocks(file)
  
  results <- data.table()
  
  for(line in taf_blocks){
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if(is.na(issued)) next
    
    # ✅ extract ALL PROBs
    matches <- str_match_all(line, "PROB(30|40).*?(\\d{4})/(\\d{4})")[[1]]
    
    if(nrow(matches) == 0) next
    
    for(j in 1:nrow(matches)){
      
      prob  <- as.numeric(matches[j,2])
      start <- matches[j,3]
      end   <- matches[j,4]
      
      start_day <- as.numeric(substr(start,1,2))
      start_hour <- as.numeric(substr(start,3,4))
      
      end_day <- as.numeric(substr(end,1,2))
      end_hour <- as.numeric(substr(end,3,4))
      
      start_time <- update(issued, day=start_day, hour=start_hour)
      end_time   <- update(issued, day=end_day, hour=end_hour)
      
      results <- rbind(results, data.table(
        Issued = issued,
        start_time = start_time,
        end_time   = end_time,
        prob       = prob
      ))
    }
  }
  
  return(results)
}

# FUNCTION 22 - EXTRACT TEMP EVENTS
extract_tempo_events <- function(file){
  
  taf <- fread(file, header=FALSE, sep="\n")
  results <- data.table()
  
  for(i in 1:nrow(taf)){
    
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if(is.na(issued)) next
    
    tws <- extract_tempo_windows(line)
    if(is.null(tws)) next
    
    for(j in 1:nrow(tws)){
      
      start_hour <- as.numeric(substr(tws$start[j],3,4))
      end_hour   <- as.numeric(substr(tws$end[j],3,4))
      
      start_time <- floor_date(issued, "day") + hours(start_hour)
      end_time   <- floor_date(issued, "day") + hours(end_hour)
      
      # ✅ FIXED TS
      has_fg <- str_detect(line,"\\bFG\\b")
      has_ts <- str_detect(line,"\\bTS\\b|TSRA|VCTS|TSGR")
      
      vis_groups <- str_extract_all(
        line,
        "(?<!\\d)\\d{4}(?!\\d)"
      )[[1]]
      
      vis_groups <- suppressWarnings(as.numeric(vis_groups))
      
      has_vis <- any(
        vis_groups <= 1000 &
          vis_groups >= 0,
        na.rm = TRUE
      )
      
      cld_vals <- str_extract_all(
        line,
        "\\b(?:BKN|OVC)\\d{3}\\b"
      )[[1]]
      
      cld_vals <- as.numeric(
        str_extract(cld_vals, "\\d{3}")
      ) * 100
      has_cld <- any(cld_vals < 1000, na.rm=TRUE)
      
      results <- rbind(results, data.table(
        Issued = issued,
        start_time = start_time,
        end_time   = end_time,
        forecast   = (has_fg | has_ts | has_vis | has_cld)
      ), fill=TRUE)
    }
  }
  
  return(results)
}

# FUNCTION 23 - Verify individual PROB events
verify_prob_events <- function(prob_events, met_dt){
  
  results <- data.table(prob = numeric(), observed = logical())
  
  for(i in 1:nrow(prob_events)){
    
    start <- prob_events$start_time[i]
    end   <- prob_events$end_time[i]
    
    obs <- met_dt[time >= start & time <= end]
    
    obs_any <- if(nrow(obs) == 0) FALSE else any(obs$event_any, na.rm=TRUE)
    
    results <- rbind(results, data.table(
      prob = prob_events$prob[i],
      observed = obs_any
    ))
  }
  
  return(results)
}

# FUNCTION 24 - BUILD TAF BLOCKS
build_taf_blocks <- function(file){
  
  taf <- fread(file, header=FALSE, sep="\n")
  
  blocks <- list()
  current_block <- ""
  
  for(i in 1:nrow(taf)){
    
    line <- taf$V1[i]
    
    # detect new TAF start (date at beginning)
    if(str_detect(line, "^\\d{8}")){
      
      if(nchar(current_block) > 0){
        blocks <- append(blocks, current_block)
      }
      
      current_block <- line
      
    } else {
      
      # continuation line
      current_block <- paste(current_block, line)
    }
  }
  
  # append final
  if(nchar(current_block) > 0){
    blocks <- append(blocks, current_block)
  }
  
  return(unlist(blocks))
}

# FUNCTION 25 - BUILD FULL DATASET
build_full_dataset <- function(taf_files, metar_files){
  
  all_data <- data.table()
  
  for(k in seq_along(taf_files)){
    
    cat("Processing:", taf_files[k], "\n")
    
    # Parse TAF (function, now with forecaster)
    taf_dt <- expand_taf_operational(taf_files[k])
    
    # Parse METAR (raw)
    met_dt_raw <- parse_metar_full(metar_files[k])
    
    # Ensure time is aligned to hour
    met_dt_raw[, time := floor_date(time, "hour")]
    
    # Collapse to one observation per hour
    met_dt <- met_dt_raw[, .(
      
      vis = median(vis,na.rm = TRUE),
      wind_spd = median(wind_spd,na.rm = TRUE),
      wind_dir = median(wind_dir, na.rm = TRUE),
      temp = median(temp,na.rm = TRUE),
      cloud_amount = median(cloud_amount,na.rm = TRUE),
      cloud_base = median(cloud_base,na.rm = TRUE),
      
      event_fg = any(event_fg),
      event_ts = any(event_ts),
      event_ra = any(event_ra),
      event_shra = any(event_shra),
      event_dz = any(event_dz),
      event_br = any(event_br),
      event_gr = any(event_gr),
      event_gs = any(event_gs)
      
      ), by = time]

    # Merge forecast + obs
    dt <- merge(
      taf_dt,
      met_dt,
      by.x = "validTime",
      by.y = "time",
      all.x = TRUE
    )
    
    print(names(dt))
    
    cat("\n")
    
    print(
      dt[
        !is.na(wind_spd_fcst)
      ][1:10]
    )
    
    # Remove rows without observations
    dt <- dt[!is.na(vis)]
    
    # Build error variables
    dt[, `:=`(
      
      fg_error =
        as.numeric(fg_fcst) -
        as.numeric(event_fg),
      
      ts_error =
        as.numeric(ts_fcst) -
        as.numeric(event_ts),
      
      wind_spd_error =
        wind_spd_fcst -
        wind_spd,
      
      wind_dir_error =
        ((wind_dir_fcst -
            wind_dir +
            180) %% 360) - 180,
      
      wind_dir_good =
        fifelse(
          
          vrb_fcst == TRUE &
            wind_spd <= 5,
          
          TRUE,
          
          abs(
            ((wind_dir_fcst -
                wind_dir +
                180) %% 360) - 180
          ) <= 20
          
        ),
      
      cloud_base_error =
        cloud_base_fcst -
        cloud_base,
      
      cloud_amount_error =
        cloud_amount_fcst -
        cloud_amount,
      
      ra_error =
        as.numeric(ra_fcst) -
        as.numeric(event_ra),
      
      shra_error =
        as.numeric(shra_fcst) -
        as.numeric(event_shra),
      
      dz_error =
        as.numeric(dz_fcst) -
        as.numeric(event_dz),
      
      br_error =
        as.numeric(br_fcst) -
        as.numeric(event_br),
      
      gr_error =
        as.numeric(gr_fcst) -
        as.numeric(event_gr),
      
      gs_error =
        as.numeric(gs_fcst) -
        as.numeric(event_gs)
    )]
    
    # Add aerodrome name
    dt[, aerodrome := toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    
    all_data <- rbind(all_data, dt, fill = TRUE)
  }
  
  return(all_data)
}

# FUNCTION 26 - CORRECT SIGNAL-BASED LEAD TIME
compute_operational_lead_times <- function(taf_dt, met_dt){
  
  results_list <- list()
  idx <- 1
  
  for(i in 1:nrow(met_dt)){
    
    t_event <- met_dt$start_time[i]
    hazard  <- met_dt$hazard[i]
    event_start <- met_dt$start_time[i]
    event_end   <- met_dt$end_time[i]
    
    # last issued TAF before event
    
    taf_before <- unique(
      taf_dt[
        Issued <= event_end &
          Issued >= event_start - hours(24),
        .(Issued)
      ]
    )
    
    if(nrow(taf_before) == 0) next
    
    setorder(taf_before, Issued)
    
    last_any <- taf_before[.N]
    
    # -----------------------
    # MATCH FORECAST SIGNAL
    # -----------------------
    
    if(hazard == "TS"){
      
      taf_signal <- taf_dt[
        ts_fcst == TRUE &
          Issued <= event_end &
          validTime >= event_start &
          validTime <= event_end
      ]
      
    } else if(hazard == "FG"){
      
      taf_signal <- taf_dt[
        fg_fcst == TRUE &
          Issued <= event_end &
          validTime >= event_start &
          validTime <= event_end
      ]
      
    } else if(hazard == "VIS"){
      
      taf_signal <- taf_dt[
        vis_fcst == TRUE &
          Issued <= event_end &
          validTime >= event_start &
          validTime <= event_end
      ]
      
    } else if(hazard == "CLD"){
      
      taf_signal <- taf_dt[
        cld1000_fcst == TRUE &
          Issued <= event_end &
          validTime >= event_start &
          validTime <= event_end
      ]
      
    }
    
    if(nrow(taf_signal) == 0){
      
      results_list[[idx]] <- data.table(
        event_time = t_event,
        hazard = hazard,
        
        first_issue = as.POSIXct(NA),
        last_issue = as.POSIXct(NA),
        last_any_issue = last_any$Issued,
        
        lead_signal_first = NA_real_,
        lead_signal_last = NA_real_,
        lead_last_any_taf =
          as.numeric(
            difftime(
              t_event,
              last_any$Issued,
              units = "hours"
            )
          ),
        
        signal_window = NA_real_,
        no_forecast = TRUE,
        late_forecast = FALSE
      )
      
      idx <- idx + 1
      next
    }
    
    setorder(taf_signal, Issued)
    
    first_signal <- taf_signal[1]
    last_signal  <- taf_signal[.N]
    
    results_list[[idx]] <- data.table(
      event_time = t_event,
      hazard = hazard,
      
      first_issue = first_signal$Issued,
      last_issue = last_signal$Issued,
      last_any_issue = last_any$Issued,
      
      lead_signal_first =
        as.numeric(
          difftime(
            t_event,
            first_signal$Issued,
            units = "hours"
          )
        ),
      
      lead_signal_last =
        as.numeric(
          difftime(
            t_event,
            last_signal$Issued,
            units = "hours"
          )
        ),
      
      lead_last_any_taf =
        as.numeric(
          difftime(
            t_event,
            last_any$Issued,
            units = "hours"
          )
        ),
      
      signal_window =
        as.numeric(
          difftime(
            last_signal$Issued,
            first_signal$Issued,
            units = "hours"
          )
        ),
      
      no_forecast = FALSE,
      late_forecast =
        first_signal$Issued > event_start &
        first_signal$Issued <= event_end
    )
    
    idx <- idx + 1
  }
  
  rbindlist(results_list, fill = TRUE)
}

# FUNCTION 27 - NEW RUNNER
run_operational_analysis <- function(taf_files, metar_files){
  
  all_results <- data.table()
  
  for(k in seq_along(taf_files)){
    
    cat("Processing:", taf_files[k], "\n")
    
    # ✅ STEP 1: expand TAF
    taf_dt <- expand_taf_operational(taf_files[k])
    
    # ✅ STEP 2: FULL METAR (no filtering)
    met_dt_full <- parse_metar_full(metar_files[k])
    
    # ✅ STEP 3: COLLAPSE INTO EVENTS
    events_dt <- collapse_events(met_dt_full)
    
    if(nrow(events_dt) == 0) next
    
    # ✅ STEP 4: COMPUTE LEAD TIMES USING EVENTS
    res <- compute_operational_lead_times(taf_dt, events_dt)
    
    if(nrow(res) == 0) next
    
    res[, aerodrome :=
          toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res, fill=TRUE)
  }
  
  return(all_results)
}

# FUNCTION 28 - COLLAPSE EVENTS TO ONE
collapse_events <- function(met_dt){
  
  # Ensure sorted
  setorder(met_dt, time)
  
  # Create hazard-specific flags
  events_list <- list()
  
  hazards <- c("event_fg", "event_ts", "event_vis_lt1000", "event_cld_lt1000")
  
  names(hazards) <- c("FG", "TS", "VIS", "CLD")
  
  for(h in names(hazards)){
    
    flag <- hazards[h]
    
    dt <- met_dt[, .(time, active = get(flag))]
    
    # Only keep rows with TRUE
    dt <- dt[active == TRUE]
    
    if(nrow(dt) == 0) next
    
    # Identify breaks in continuity (gap > 1 hour)
    dt[, gap := as.numeric(difftime(time, shift(time), units="hours"))]
    dt[is.na(gap), gap := 0]
    
    # New event if gap > 1 hour
    dt[, new_event := (gap > 1)]
    
    dt[, event_id := cumsum(new_event) + 1]
    
    # Collapse to one row per event
    collapsed <- dt[, .(
      start_time = min(time),
      end_time   = max(time),
      duration_hours =
        as.numeric(
          difftime(
            max(time),
            min(time),
            units = "hours"
          )
        ) + 1
    ), by = event_id]
    
    collapsed[, hazard := h]
    
    events_list[[h]] <- collapsed
  }
  
  events <- rbindlist(events_list, fill=TRUE)
  
  return(events)
}

# FUNCTION 29
extract_taf_blocks <- function(line){
  
  starts <- gregexpr(
    "FM\\d{6}|BECMG\\s+\\d{4}/\\d{4}|TEMPO\\s+\\d{4}/\\d{4}|INTER\\s+\\d{4}/\\d{4}|PROB30\\s+\\d{4}/\\d{4}|PROB40\\s+\\d{4}/\\d{4}",
    line,
    perl = TRUE
  )[[1]]
  
  if(starts[1] == -1){
    return(data.table(
      type = "BASE",
      text = line
    ))
  }
  
  starts <- c(1, starts)
  ends <- c(starts[-1] - 1, nchar(line))
  
  data.table(
    text = substring(line, starts, ends)
  )[
    ,
    type := fifelse(
      grepl("^FM", text), "FM",
      fifelse(
        grepl("^BECMG", text), "BECMG",
        fifelse(
          grepl("^TEMPO", text), "TEMPO",
          fifelse(
            grepl("^INTER", text), "INTER",
            fifelse(
              grepl("^PROB30", text), "PROB30",
              fifelse(
                grepl("^PROB40", text),
                "PROB40",
                "BASE"
              )
            )
          )
        )
      )
    )
  ]
}

# FUNCTION 30
extract_block_hazards <- function(txt){
  
  # remove remarks
  txt <- strsplit(txt, "\\bRMK\\b")[[1]][1]
  
  ts_flag <- str_detect(
    txt,
    "\\bTS\\b|TSRA|VCTS|TSGR"
  )
  
  fg_flag <- str_detect(
    txt,
    "FG|BCFG|MIFG|PRFG"
  )
  
  vis_groups <- str_match_all(
    txt,
    "\\b(\\d{4})\\s+(?:FG|BR|HZ|FU|DZ|RA|SHRA|TSRA|VCTS|FEW|SCT|BKN|OVC)"
  )[[1]]
  
  vis_vals <- if(nrow(vis_groups) > 0){
    as.numeric(vis_groups[,2])
  } else {
    numeric(0)
  }
  
  vis_flag <- any(vis_vals < 5000, na.rm = TRUE)
  
  # -----------------------------
  # Cloud information
  # -----------------------------
  
  cloud_info <- extract_cloud_layers(txt)
  
  cld_flag <- cloud_info$low_cloud
  
  # -----------------------------
  # Wind
  # -----------------------------
  
  wind_match <- str_match(
    txt,
    "\\b(\\d{3}|VRB)(\\d{2})KT\\b"
  )
  
  wind_dir_txt <- wind_match[2]
  
  vrb_fcst <- identical(
    wind_dir_txt,
    "VRB"
  )
  
  wind_dir <- suppressWarnings(
    as.numeric(wind_dir_txt)
  )
  
  wind_spd <- suppressWarnings(
    as.numeric(wind_match[3])
  )
  
  if(length(wind_dir)==0) wind_dir <- NA_real_
  
  if(length(wind_spd)==0) wind_spd <- NA_real_
  
  # -----------------------------
  # Weather phenomena
  # -----------------------------
  
  wx_info <- extract_weather_phenomena(
    txt
  )
  
  
  list(
    fg = fg_flag,
    ts = ts_flag,
    vis = vis_flag,
    cld = cld_flag,
    wind_dir = wind_dir,
    wind_spd = wind_spd,
    vrb_fcst = vrb_fcst,
    cloud_amount = cloud_info$cloud_amount,
    cloud_base = cloud_info$cloud_base,
    ra = wx_info$ra,
    shra = wx_info$shra,
    dz = wx_info$dz,
    br = wx_info$br,
    gr = wx_info$gr,
    gs = wx_info$gs
  )
}

# FUNCTION 31
decode_ddhh <- function(issued, ddhh){
  
  day  <- as.numeric(substr(ddhh,1,2))
  hour <- as.numeric(substr(ddhh,3,4))
  
  yr  <- year(issued)
  mon <- month(issued)
  
  candidate <- as.POSIXct(
    sprintf(
      "%04d-%02d-%02d %02d:00:00",
      yr,
      mon,
      day,
      hour
    ),
    tz = "UTC"
  )
  
  candidate
}

# FUNCTION 32
extract_taf_temperature <- function(
    txt,
    valid_start,
    valid_end
){
  
  temp_match <- str_match(
    txt,
    "RMK\\s+T\\s+([0-9\\s\\-]+)"
  )
  
  if(is.na(temp_match[2])){
    return(NULL)
  }
  
  temps <- suppressWarnings(
    as.numeric(
      strsplit(
        trimws(temp_match[2]),
        "\\s+"
      )[[1]]
    )
  )
  
  if(length(temps) == 0){
    return(NULL)
  }
  
  temp_times <- valid_start + hours(
    seq_len(length(temps)) * 3
  )
  
  data.table(
    validTime = temp_times,
    temp_fcst = temps
  )
}

# FUNCTION 33 - TAF OPERATIONAL MASTER FUNCTION
expand_taf_operational <- function(file){
  
  taf_blocks <- build_taf_blocks(file)
  
  out <- data.table()
  
  for(line in taf_blocks){
    
    date_str <- str_extract(line, "^\\d{8}")
    
    issue_time <- str_extract(
      line,
      "(?<=\\s)\\d{6}(?=\\s)"
    )
    
    issued <- suppressWarnings(
      ymd_hms(
        paste(date_str, issue_time),
        tz = "UTC"
      )
    )
    
    if(is.na(issued)) next
    
    # -----------------------------
    # FORECASTER
    # -----------------------------
    forecaster <- str_match(
      line,
      "^\\d{8}\\s+\\d{6}\\s+(\\S+)"
    )[,2]
    
    if(is.na(forecaster))
      forecaster <- "UNKNOWN"
    
    # -----------------------------
    # VALIDITY PERIOD
    # -----------------------------
    
    valid_period <- str_match(
      line,
      "\\b(\\d{4})\\s+(\\d{4})\\b"
    )
    
    if(any(is.na(valid_period))) next
    
    valid_start <- decode_ddhh(
      issued,
      valid_period[2]
    )
    
    valid_end <- decode_ddhh(
      issued,
      valid_period[3]
    )
    
    if(valid_end < valid_start){
      valid_end <- valid_end %m+% months(1)
    }
    
    dt <- data.table(
      Issued = issued,
      
      validTime = seq(
        valid_start,
        valid_end,
        by = "1 hour"
      ),
      
      forecaster = forecaster,
      
      forecast_type = "BASE",
      
      cloud_amount_fcst = NA_real_,
      cloud_base_fcst = NA_real_,
      
      ra_fcst = FALSE,
      shra_fcst = FALSE,
      dz_fcst = FALSE,
      br_fcst = FALSE,
      gr_fcst = FALSE,
      gs_fcst = FALSE,
      
      fg_fcst = FALSE,
      ts_fcst = FALSE,
      vis_fcst = FALSE,
      cld1000_fcst = FALSE,
      
      wind_dir_fcst = NA_real_,
      wind_spd_fcst = NA_real_,
      vrb_fcst = FALSE
    )
    
    # -----------------------------
    # BASE FORECAST
    # -----------------------------
    
    base_text <- strsplit(
      line,
      "FM\\d{6}|BECMG\\s+\\d{4}/\\d{4}|TEMPO\\s+\\d{4}/\\d{4}|INTER\\s+\\d{4}/\\d{4}|PROB30\\s+\\d{4}/\\d{4}|PROB40\\s+\\d{4}/\\d{4}"
    )[[1]][1]
    
    base_hz <- extract_block_hazards(
      base_text
    )
    
    dt[, `:=`(
      fg_fcst = base_hz$fg,
      ts_fcst = base_hz$ts,
      vis_fcst = base_hz$vis,
      cld1000_fcst = base_hz$cld,
      
      wind_dir_fcst = base_hz$wind_dir,
      wind_spd_fcst = base_hz$wind_spd,
      vrb_fcst = base_hz$vrb_fcst,
      
      cloud_amount_fcst = base_hz$cloud_amount,
      cloud_base_fcst = base_hz$cloud_base,
      ra_fcst = base_hz$ra,
      shra_fcst = base_hz$shra,
      dz_fcst = base_hz$dz,
      br_fcst = base_hz$br,
      gr_fcst = base_hz$gr,
      gs_fcst = base_hz$gs
    )]
    
    # ===================================================
    # FM GROUPS
    # ===================================================
    
    fm_matches <- str_match_all(
      line,
      "FM(\\d{6})"
    )[[1]]
    
    if(nrow(fm_matches) > 0){
      
      fm_pos <- gregexpr(
        "FM\\d{6}",
        line
      )[[1]]
      
      for(j in seq_len(nrow(fm_matches))){
        
        fm_code <- fm_matches[j,2]
        
        start_time <- decode_ddhh(
          issued,
          substr(fm_code,1,4)
        )
        
        text_start <- fm_pos[j]
        
        if(j < nrow(fm_matches)){
          text_end <- fm_pos[j + 1] - 1
          
          next_start <- decode_ddhh(
            issued,
            substr(
              fm_matches[j + 1,2],
              1,
              4
            )
          )
          
          end_time <- next_start - hours(1)
          
        } else {
          
          text_end <- nchar(line)
          end_time <- valid_end
        }
        
        fm_text <- substr(
          line,
          text_start,
          text_end
        )
        
        hz <- extract_block_hazards(
          fm_text
        )
        
        dt[
          validTime >= start_time &
            validTime <= end_time,
          `:=`(
            forecast_type = "FM",
            
            fg_fcst = hz$fg,
            ts_fcst = hz$ts,
            vis_fcst = hz$vis,
            cld1000_fcst = hz$cld,
            
            wind_dir_fcst = hz$wind_dir,
            wind_spd_fcst = hz$wind_spd,
            vrb_fcst = hz$vrb_fcst,
            
            cloud_amount_fcst = hz$cloud_amount,
            cloud_base_fcst = hz$cloud_base,
            ra_fcst = hz$ra,
            shra_fcst = hz$shra,
            dz_fcst = hz$dz,
            br_fcst = hz$br,
            gr_fcst = hz$gr,
            gs_fcst = hz$gs
          )
        ]
      }
    }
    
    # ===================================================
    # BECMG
    # ===================================================
    
    becmg_matches <- str_match_all(
      line,
      "BECMG\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if(nrow(becmg_matches) > 0){
      
      for(j in seq_len(nrow(becmg_matches))){
        
        end_time <- decode_ddhh(
          issued,
          becmg_matches[j,3]
        )
        
        hz <- extract_block_hazards(
          becmg_matches[j,4]
        )
        
        dt[
          validTime >= end_time,
          `:=`(
            forecast_type = "BECMG",
            
            fg_fcst = hz$fg,
            ts_fcst = hz$ts,
            vis_fcst = hz$vis,
            cld1000_fcst = hz$cld,
            
            wind_dir_fcst = hz$wind_dir,
            wind_spd_fcst = hz$wind_spd,
            vrb_fcst = hz$vrb_fcst,
            
            cloud_amount_fcst = hz$cloud_amount,
            cloud_base_fcst = hz$cloud_base,
            ra_fcst = hz$ra,
            shra_fcst = hz$shra,
            dz_fcst = hz$dz,
            br_fcst = hz$br,
            gr_fcst = hz$gr,
            gs_fcst = hz$gs
          )
        ]
      }
    }
    
    # ===================================================
    # TEMPO
    # ===================================================
    
    tempo_matches <- str_match_all(
      line,
      "TEMPO\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if(nrow(tempo_matches) > 0){
      
      for(j in seq_len(nrow(tempo_matches))){
        
        start_time <- decode_ddhh(
          issued,
          tempo_matches[j,2]
        )
        
        end_time <- decode_ddhh(
          issued,
          tempo_matches[j,3]
        )
        
        hz <- extract_block_hazards(
          tempo_matches[j,4]
        )
        
        dt[
          validTime >= start_time &
            validTime <= end_time,
          `:=`(
            forecast_type = "TEMPO",
            
            fg_fcst = hz$fg,
            ts_fcst = hz$ts,
            vis_fcst = hz$vis,
            cld1000_fcst = hz$cld,
            
            wind_dir_fcst = hz$wind_dir,
            wind_spd_fcst = hz$wind_spd,
            vrb_fcst = hz$vrb_fcst,
            
            cloud_amount_fcst = hz$cloud_amount,
            cloud_base_fcst = hz$cloud_base,
            ra_fcst = hz$ra,
            shra_fcst = hz$shra,
            dz_fcst = hz$dz,
            br_fcst = hz$br,
            gr_fcst = hz$gr,
            gs_fcst = hz$gs
          )
        ]
      }
    }
    
    # ===================================================
    # INTER
    # ===================================================
    
    inter_matches <- str_match_all(
      line,
      "INTER\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if(nrow(inter_matches) > 0){
      
      for(j in seq_len(nrow(inter_matches))){
        
        start_time <- decode_ddhh(
          issued,
          inter_matches[j,2]
        )
        
        end_time <- decode_ddhh(
          issued,
          inter_matches[j,3]
        )
        
        hz <- extract_block_hazards(
          inter_matches[j,4]
        )
        
        dt[
          validTime >= start_time &
            validTime <= end_time,
          `:=`(
            forecast_type = "INTER",
            
            fg_fcst = hz$fg,
            ts_fcst = hz$ts,
            vis_fcst = hz$vis,
            cld1000_fcst = hz$cld,
            
            wind_dir_fcst = hz$wind_dir,
            wind_spd_fcst = hz$wind_spd,
            vrb_fcst = hz$vrb_fcst,
            
            cloud_amount_fcst = hz$cloud_amount,
            cloud_base_fcst = hz$cloud_base,
            ra_fcst = hz$ra,
            shra_fcst = hz$shra,
            dz_fcst = hz$dz,
            br_fcst = hz$br,
            gr_fcst = hz$gr,
            gs_fcst = hz$gs
          )
        ]
      }
    }
    
    # ===================================================
    # PROB30 / PROB40
    # ===================================================
    
    prob_matches <- str_match_all(
      line,
      "PROB(30|40)\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if(nrow(prob_matches) > 0){
      
      for(j in seq_len(nrow(prob_matches))){
        
        start_time <- decode_ddhh(
          issued,
          prob_matches[j,3]
        )
        
        end_time <- decode_ddhh(
          issued,
          prob_matches[j,4]
        )
        
        hz <- extract_block_hazards(
          prob_matches[j,5]
        )
        
        dt[
          validTime >= start_time &
            validTime <= end_time,
          `:=`(
            forecast_type =
              paste0(
                "PROB",
                prob_matches[j,2]
              ),
            
            fg_fcst = hz$fg,
            ts_fcst = hz$ts,
            vis_fcst = hz$vis,
            cld1000_fcst = hz$cld,
            
            wind_dir_fcst = hz$wind_dir,
            wind_spd_fcst = hz$wind_spd,
            vrb_fcst = hz$vrb_fcst,
            
            cloud_amount_fcst = hz$cloud_amount,
            cloud_base_fcst = hz$cloud_base,
            ra_fcst = hz$ra,
            shra_fcst = hz$shra,
            dz_fcst = hz$dz,
            br_fcst = hz$br,
            gr_fcst = hz$gr,
            gs_fcst = hz$gs
          )
        ]
      }
    }
    
    out <- rbind(
      out,
      dt,
      fill = TRUE
    )
  }
  
  out
}

# FUNCTION 34 - EXTRACT TEMP AND VALIDATE NOISE MASTER FUNCTION
extract_taf_temperatures <- function(file){
  
  taf_blocks <- build_taf_blocks(file)
  
  out <- data.table()
  
  for(line in taf_blocks){
    
    date_str <- str_extract(
      line,
      "^\\d{8}"
    )
    
    issue_time <- str_extract(
      line,
      "(?<=\\s)\\d{6}(?=\\s)"
    )
    
    issued <- suppressWarnings(
      ymd_hms(
        paste(date_str, issue_time),
        tz = "UTC"
      )
    )
    
    if(is.na(issued))
      next
    
    forecaster <- str_match(
      line,
      "^\\d{8}\\s+\\d{6}\\s+(\\S+)"
    )[,2]
    
    valid_period <- str_match(
      line,
      "\\b(\\d{4})\\s+(\\d{4})\\b"
    )
    
    if(any(is.na(valid_period)))
      next
    
    valid_start <- decode_ddhh(
      issued,
      valid_period[2]
    )
    
    # -------------------------
    # RMK T extraction
    # -------------------------
    
    t_match <- str_match(
      line,
      "RMK\\s+T\\s+([0-9\\-]+)\\s+([0-9\\-]+)\\s+([0-9\\-]+)\\s+([0-9\\-]+)"
    )
    
    if(all(is.na(t_match)))
      next
    
    temps <- as.numeric(
      t_match[2:5]
    )
    
    temp_times <- valid_start +
      hours(c(3,6,9,12))
    
    dt <- data.table(
      Issued = issued,
      forecaster = forecaster,
      validTime = temp_times,
      temp_fcst = temps
    )
    
    out <- rbind(
      out,
      dt,
      fill = TRUE
    )
  }
  
  out
}

# FUNCTION 35 - BUILD TEMP DATA MASTER FUNCTION
build_temperature_dataset <- function(
    taf_files,
    metar_files
){
  
  all_data <- data.table()
  
  for(k in seq_along(taf_files)){
    
    taf_temp <- extract_taf_temperatures(
      taf_files[k]
    )
    
    met <- parse_metar_full(
      metar_files[k]
    )
    
    met <- met[
      ,
      .(
        temp = median(
          temp,
          na.rm=TRUE
        )
      ),
      by=time
    ]
    
    dt <- merge(
      taf_temp,
      met,
      by.x="validTime",
      by.y="time",
      all.x=TRUE
    )
    
    dt[, temp_error :=
         temp_fcst - temp]
    
    dt[, aerodrome :=
         toupper(
           gsub(
             "TAF_|_2025\\.txt",
             "",
             taf_files[k]
           )
         )]
    
    all_data <- rbind(
      all_data,
      dt,
      fill=TRUE
    )
  }
  
  all_data
}

# FUNCTION 36 - CLOUD LAYER EXTRACTION
extract_cloud_layers <- function(txt){
  
  cloud_match <- str_match_all(
    txt,
    "\\b(FEW|SCT|BKN|OVC)(\\d{3})\\b"
  )[[1]]
  
  if(nrow(cloud_match) == 0){
    return(
      list(
        cloud_amount = NA_real_,
        cloud_base = NA_real_,
        low_cloud = FALSE
      )
    )
  }
  
  amount_lookup <- c(
    FEW = 1,
    SCT = 2,
    BKN = 3,
    OVC = 4
  )
  
  cloud_amount <- max(
    amount_lookup[
      cloud_match[,2]
    ],
    na.rm = TRUE
  )
  
  cloud_base <- min(
    as.numeric(
      cloud_match[,3]
    ) * 100,
    na.rm = TRUE
  )
  
  low_cloud <- any(
    cloud_match[,2] %in%
      c("BKN","OVC") &
      (as.numeric(
        cloud_match[,3]
      ) * 100 < 1000)
  )
  
  list(
    cloud_amount = cloud_amount,
    cloud_base = cloud_base,
    low_cloud = low_cloud
  )
}

# FUNCTION 37 - WEATHER PHENOMENA EXTRACTION
extract_weather_phenomena <- function(txt){
  
  list(
    
    ra =
      grepl(
        "\\bRA\\b",
        txt
      ),
    
    shra =
      grepl(
        "SHRA",
        txt
      ),
    
    dz =
      grepl(
        "\\bDZ\\b",
        txt
      ),
    
    br =
      grepl(
        "\\bBR\\b",
        txt
      ),
    
    fg =
      grepl(
        "\\bFG\\b",
        txt
      ),
    
    ts =
      grepl(
        "TS",
        txt
      ),
    
    gr =
      grepl(
        "\\bGR\\b",
        txt
      ),
    
    gs =
      grepl(
        "\\bGS\\b",
        txt
      )
    
  )
}

# FUNCTION 38 - BUILD COMPOSITE NOISE INDEX
build_composite_noise_index <- function(
    noise_data,
    temp_data
){
  
  temp_metric <- temp_data[
    ,
    .(
      temp_mae =
        mean(
          abs(temp_error),
          na.rm = TRUE
        )
    ),
    by = ForecasterID
  ]
  
  wind_metric <- noise_data[
    ,
    .(
      
      wind_spd_mae =
        mean(
          abs(wind_spd_error),
          na.rm = TRUE
        ),
      
      wind_dir_mae =
        mean(
          abs(wind_dir_error),
          na.rm = TRUE
        )
      
    ),
    by = ForecasterID
  ]
  
  cloud_metric <- noise_data[
    ,
    .(
      
      cloud_base_mae =
        mean(
          abs(
            cloud_base_fcst -
              cloud_base
          ),
          na.rm = TRUE
        ),
      
      cloud_amount_mae =
        mean(
          abs(
            cloud_amount_fcst -
              cloud_amount
          ),
          na.rm = TRUE
        )
      
    ),
    by = ForecasterID
  ]
  
  hazard_metric <- noise_data[
    ,
    .(
      
      fg_rate =
        mean(fg_fcst),
      
      ts_rate =
        mean(ts_fcst),
      
      vis_rate =
        mean(vis_fcst),
      
      cld_rate =
        mean(cld1000_fcst)
      
    ),
    by = ForecasterID
  ]
  
  structure_metric <- structure_entropy
  
  out <- Reduce(
    
    function(x,y)
      merge(
        x,
        y,
        by="ForecasterID",
        all=TRUE
      ),
    
    list(
      
      temp_metric,
      
      wind_metric,
      
      cloud_metric,
      
      hazard_metric,
      
      structure_metric
      
    )
    
  )
  
  metric_cols <-
    
    names(out)[
      names(out) != "ForecasterID"
    ]
  
  for(col in metric_cols){
    
    out[
      ,
      paste0(
        col,
        "_z"
      ) := scale(
        get(col)
      )
    ]
    
  }
  
  z_cols <-
    
    grep(
      "_z$",
      names(out),
      value = TRUE
    )
  
  out[
    ,
    CompositeNoiseIndex :=
      rowMeans(
        .SD,
        na.rm = TRUE
      ),
    .SDcols = z_cols
  ]
  
  out[
    order(
      -CompositeNoiseIndex
    )
  ]
}


#######################################
# Load data
### READ TAF DATA
### AVS1 temp and wind
### Data available from the author on request
taf_files<-list.files(pattern = "TAF")
metar_files<-list.files(pattern = "METAR")
taf_files<-sort(taf_files)
metar_files<-sort(metar_files)

# test
# test <- expand_taf_with_amd(taf_files[1])
# unique(test$forecaster)

# Run BUILD DATASET function
noise_data <- build_full_dataset(taf_files, metar_files)
temp_data <- build_temperature_dataset(taf_files,metar_files)

# Immediate check
summary(temp_data$temp_fcst)
summary(temp_data$temp)
summary(temp_data$temp_error)

# Anonymise Forecasters
set.seed(2025)
forecaster_lookup <- data.table(
  forecaster =sort(unique(noise_data$forecaster))
  )

forecaster_lookup[,
  ForecasterID := sample(sprintf("F%03d",seq_len(.N)))
  ]

noise_data <- merge(noise_data,forecaster_lookup,by="forecaster",all.x=TRUE)
temp_data <- merge(temp_data,forecaster_lookup,by="forecaster",all.x=TRUE)

#############################################
### NOISE AUDIT ACTIONS
#############################################
# Part 1 — Hazard Propensity Noise
noise_data[, n_forecaster := .N, by = forecaster]
noise_data_filt <- noise_data[n_forecaster >= 30]

hazard_noise <- noise_data_filt[
  ,
  .(
    fg_rate = mean(fg_fcst),
    ts_rate = mean(ts_fcst),
    vis_rate = mean(vis_fcst),
    cld_rate = mean(cld1000_fcst)
  ),
  by = ForecasterID
]

hazard_noise

hazard_noise_summary <- hazard_noise[
  ,
  .(
    mean_fg = mean(fg_rate),
    sd_fg = sd(fg_rate),
    mean_ts = mean(ts_rate),
    sd_ts = sd(ts_rate),
    mean_vis = mean(vis_rate),
    sd_vis = sd(vis_rate),
    mean_cld = mean(cld_rate),
    sd_cld = sd(cld_rate)
  )
]

hazard_noise_summary

# Plot hazard Noise
ggplot(melt(hazard_noise,id.vars="ForecasterID"),aes(value)) +
  geom_histogram() +
  facet_wrap(~variable)

hazard_long <- melt(
  hazard_noise,
  id.vars = "ForecasterID",
  variable.name = "Hazard",
  value.name = "ForecastRate"
)

hazard_long[
  ,
  Hazard := factor(
    Hazard,
    levels = c(
      "fg_rate",
      "ts_rate",
      "vis_rate",
      "cld_rate"
    ),
    labels = c(
      "Fog (FG)",
      "Thunderstorms (TS)",
      "Low Visibility (VIS)",
      "Low Cloud (CLD)"
    )
  )
]

ggplot(
  hazard_long,
  aes(
    x = Hazard,
    y = ForecastRate
  )
) +
  geom_violin(
    fill = "lightblue",
    colour = "black",
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white"
  ) +
  labs(
    title = "Hazard Propensity Noise",
    subtitle = "Distribution of hazard forecast frequencies across forecasters",
    x = NULL,
    y = "Forecast Frequency",
    caption = "Each point represents a forecaster's forecast propensity."
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold")
  )

#################################################
# Part 2 — Airport-Adjusted Level Noise
hazard_airport <- noise_data_filt[
  ,
  .(
    n = .N,
    ts_rate = mean(ts_fcst),
    fg_rate = mean(fg_fcst),
    vis_rate = mean(vis_fcst),
    cld_rate = mean(cld1000_fcst)
  ),
  by = .(ForecasterID,aerodrome)
]

hazard_airport <- hazard_airport[
  n >= 30
]

hazard_airport

hazard_airport[
  ,
  ts_anomaly := ts_rate - mean(ts_rate),
  by = aerodrome
]

hazard_airport

hazard_airport[
  ,
  fg_anomaly := fg_rate - mean(fg_rate),
  by = aerodrome
]

hazard_airport

hazard_airport[
  ,
  vis_anomaly := vis_rate - mean(vis_rate),
  by = aerodrome
]

hazard_airport

hazard_airport[
  ,
  cld_anomaly := cld_rate - mean(cld_rate),
  by = aerodrome
]

hazard_airport

hazard_long <- melt(hazard_airport, 
                    measure.vars = c("ts_anomaly","fg_anomaly","vis_anomaly","cld_anomaly")
                    )

# Plot
ggplot(hazard_long,aes(variable,value)) +
  geom_boxplot()

ni_labels <- data.table(
  Hazard = c(
    "Fog",
    "Thunderstorms",
    "Low Visibility",
    "Low Cloud"
  ),
  label = c(
    "NI = 0.103",
    "NI = 0.037",
    "NI = 0.096",
    "NI = 0.151"
  ),
  y = c(0.30,0.30,0.30,0.30)
)

hazard_long <- melt(
  hazard_airport,
  measure.vars = c(
    "fg_anomaly",
    "ts_anomaly",
    "vis_anomaly",
    "cld_anomaly"
  ),
  variable.name = "Hazard",
  value.name = "Anomaly"
)

hazard_long[
  ,
  Hazard := factor(
    Hazard,
    levels = c(
      "fg_anomaly",
      "ts_anomaly",
      "vis_anomaly",
      "cld_anomaly"
    ),
    labels = c(
      "Fog",
      "Thunderstorms",
      "Low Visibility",
      "Low Cloud"
    )
  )
]

ggplot(hazard_long,aes(x = Hazard,y = Anomaly,fill = Hazard)) +
  geom_violin(alpha = 0.5,colour = NA) +
  geom_boxplot(width = 0.12,outlier.size = 0.8,fill = "white") +
  geom_hline(yintercept = 0,linetype = 2,colour = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Airport-Adjusted Forecaster Noise",
       subtitle = "Hazard forecast anomalies relative to airport peer averages",
       x = NULL,
       y = "Forecast Propensity Anomaly",
    caption =
      paste(
        "Positive values = more frequent forecasting than peers at same location;",
        "negative values indicate less.")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold",size = 16),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(face = "bold")) +
  geom_text(data = ni_labels,aes(x = Hazard,y = y,label = label),
            inherit.aes = FALSE,
            fontface = "bold",
            size = 4) + ylim(c(-0.25,0.4))

############################################
# Table 2
level_noise_summary <- hazard_airport[
  ,
  .(
    ts_noise = sd(ts_anomaly),
    fg_noise = sd(fg_anomaly),
    vis_noise = sd(vis_anomaly),
    cld_noise = sd(cld_anomaly)
  )
]

level_noise_summary

# convert to noise index:
hazard_NI <- data.table(
  hazard =
    c("TS","FG","VIS","CLD"),
  tau =
    c(
      level_noise_summary$ts_noise,
      level_noise_summary$fg_noise,
      level_noise_summary$vis_noise,
      level_noise_summary$cld_noise
    )
)

hazard_NI[
  ,
  NI :=
    (2/sqrt(pi))*tau
]

hazard_NI

##########################################
## Figure 3
structure_noise <- noise_data_filt[
  !forecast_type %in% c("BASE","FM"),
  .N,
  by = .(ForecasterID,forecast_type)
]

structure_noise

structure_noise[, pct := N / sum(N), by = ForecasterID]

structure_noise

ggplot(structure_noise,aes(ForecasterID,pct,fill=forecast_type)) +
  geom_col()


#############################################
## Figure 4 - Structure Entropy
structure_entropy <- structure_noise[
  ,
  .(
    entropy =
      -sum(
        pct * log(pct),
        na.rm = TRUE
      )
  ),
  by = ForecasterID
]

structure_entropy[
  order(-entropy)
]

ggplot(structure_entropy[order(entropy)],
  aes(reorder(ForecasterID,entropy),entropy)) +
  geom_col()

top_entropy <- structure_entropy[
  order(-entropy)
]

# Plot
ggplot(top_entropy,aes(reorder(ForecasterID,entropy),entropy)) +
  geom_col(fill = "darkgreen") +
  coord_flip()

top_entropy_nonzero <- top_entropy[entropy > 0]

entropy_plot <- rbind(
  top_entropy_nonzero[1:15],
  top_entropy_nonzero[(.N-14):.N]
)

ggplot(entropy_plot,aes(reorder(ForecasterID, entropy),entropy)) +
  geom_col(aes(fill = entropy)) +
  coord_flip() +
  scale_fill_viridis_c(guide = "none") +
  labs(title = "Forecast Communication Diversity",
       subtitle = "Top and bottom 15 forecasters ranked by forecast structure entropy",
       x = "Forecaster ID",
       y = "Communication Entropy") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold",size = 16),
        axis.title = element_text(face = "bold"))

#######################
# Figure 5 - Structure Noise Summary
structure_summary <- structure_noise[
  ,
  .(
    mean_pct = mean(pct),
    sd_pct = sd(pct)
  ),
  by = forecast_type
]

structure_summary[order(-sd_pct)]

ggplot(structure_summary, aes(reorder(forecast_type,sd_pct),sd_pct)) +
  geom_col()

# Plot
ggplot(structure_summary,aes(reorder(forecast_type,sd_pct),sd_pct)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Forecast Structure Noise",
       subtitle ="Inter-forecaster variability in forecast structure selection",
       x = NULL, 
       y = "Standard Deviation of Usage Rate",
       caption = "Higher values indicate greater variability in how forecasters communicate forecast evolution and uncertainty."
  ) +
  theme_minimal(base_size = 13)


################################################
# Part 4 — Temperature Noise
# Build Accuracy Table
temp_audit <- temp_data[
  !is.na(temp_error),
  
  .(
    n=.N,
    accuracy =
      mean(
        abs(temp_error)<=1
      ),
    mae =
      mean(
        abs(temp_error)
      ),
    bias =
      mean(
        temp_error
      )
  ),
  by=ForecasterID
]

# Figure 6
ggplot(temp_audit,aes(reorder(ForecasterID,accuracy),accuracy)) +
  geom_col()

# Plot
ggplot(temp_audit,aes(reorder(ForecasterID, accuracy),accuracy)) +
  geom_col(fill = "#2C7FB8") +
  coord_flip() +
  geom_hline(
    yintercept = mean(temp_audit$accuracy),
    colour = "red",
    linetype = 2,
    linewidth = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits = c(0, 1)) +
  labs(
    title =
      "Temperature Forecast Accuracy by Forecaster",
    subtitle =
      expression(
        paste(
          "Accuracy defined as |forecast error| ",
          le,
          " 1" * degree * "C"
        )
      ),
    x = "Forecaster ID",
    y = "Accuracy",
    caption = "Bars show the proportion of temperature forecasts within ±1°C of the corresponding observation. The red dashed line shows the mean accuracy across all forecasters."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold",size = 16),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_blank()
  )

ggplot(temp_audit, aes(x = accuracy, y = reorder(ForecasterID, accuracy))) +
  geom_point(size = 2.5, colour = "#2C7FB8") +
  geom_vline(
    xintercept = mean(temp_audit$accuracy),
    colour = "red",
    linetype = 2
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title =
      "Temperature Forecast Accuracy by Forecaster",
    subtitle =
      "Forecasters ranked by proportion of forecasts within ±1°C",
    x = "Accuracy",
    y = "Forecaster ID",
    caption = "Red dashed line indicates the mean forecaster accuracy."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title =
          element_text(face = "bold", size = 16),
        axis.text.y =
          element_text(size = 6))

######################################
# Part 5 — Wind Speed Noise
wspd_audit <- noise_data[
  !is.na(wind_spd_error),
  
  .(
    n = .N,
    accuracy =
      mean(
        abs(wind_spd_error) <= 5
      ),
    mae =
      mean(
        abs(wind_spd_error)
      ),
    bias =
      mean(
        wind_spd_error
      )
  ),
  by = ForecasterID
]

# Figure
ggplot(wspd_audit,aes(reorder(ForecasterID,accuracy),accuracy)) +
  geom_col()

### Dual Wind Temp Plots
temp_plot <- temp_audit[
  ,
  .(
    ForecasterID,
    accuracy,
    Variable = "Temperature"
  )
]
wspd_plot <- wspd_audit[
  ,
  .(
    ForecasterID,
    accuracy,
    Variable = "Wind Speed"
  )
]
acc_plot <- rbind(
  temp_plot,
  wspd_plot
)
#Plot
ggplot(
  acc_plot,
  aes(
    x = accuracy,
    y = reorder(
      ForecasterID,
      accuracy
    ),
    colour = Variable
  )
) +
  geom_point(
    size = 2.7,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free_y"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    title =
      "Forecast Accuracy Variability by Forecaster",
    subtitle =
      paste(
        "Temperature accuracy defined as ±1°C;",
        "Wind speed accuracy defined as ±5 kt"
      ),
    x = "Forecast Accuracy",
    y = "Forecaster ID",
    caption =
      paste(
        "Forecasters are ranked independently within each forecast variable.",
        "Greater spread indicates higher forecast noise."
      )
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    plot.subtitle =
      element_text(
        size = 11
      ),
    axis.title =
      element_text(
        face = "bold"
      ),
    axis.text.y =
      element_text(
        size = 6
      ),
    legend.position = "none"
  )

# Plot
acc_plot <- rbind(
  temp_audit[
    ,
    .(
      ForecasterID,
      accuracy,
      Variable = "Temperature"
    )
  ],
  wspd_audit[
    ,
    .(
      ForecasterID,
      accuracy,
      Variable = "Wind Speed"
    )
  ]
)

ggplot(
  acc_plot,
  aes(
    Variable,
    accuracy,
    fill = Variable
  )
) +
  geom_violin(
    alpha = 0.4
  ) +
  
  geom_boxplot(
    width = 0.12,
    fill = "white"
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    title =
      "Forecast Accuracy Noise",
    
    subtitle =
      "Distribution of forecaster accuracies for temperature and wind speed forecasts",
    
    x = NULL,
    
    y =
      "Forecast Accuracy",
    
    caption =
      "Temperature accuracy forecasts within ±1°C; wind speed accuracy forecasts within ±5 kt."
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    
    axis.title =
      element_text(
        face = "bold"
      ),
    
    legend.position = "none"
  )

##################################################
### Part 6 — Wind Direction Noise
wdir_audit <- noise_data[
  !is.na(wind_dir_good),
  .(
    n = .N,
    accuracy =
      mean(
        wind_dir_good
      )
  ),
  by = ForecasterID
]

wdir_audit

# Figure 8
ggplot(wdir_audit,aes(reorder(ForecasterID,accuracy),accuracy)) +
  geom_col()

### COMNINED VOLIN PLOT
acc_plot <- rbind(temp_audit[, .(ForecasterID, accuracy, Variable = "Temperature")], 
                  wspd_audit[, .(ForecasterID, accuracy, Variable = "Wind Speed")], 
                  wdir_audit[, .(ForecasterID, accuracy, Variable = "Wind Direction")])

## PLOT
ggplot(
  acc_plot,
  aes(
    Variable,
    accuracy,
    fill = Variable
  )
) +
  
  geom_violin(
    alpha = 0.4,
    colour = NA
  ) +
  
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.size = 0.8
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "red",
    colour = "black"
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0,1)
  ) +
  
  labs(
    title =
      "Continuous Forecast Accuracy Noise",
    
    subtitle =
      paste(
        "Distribution of forecaster accuracies for temperature, wind speed",
        "and wind direction forecasts"
      ),
    
    x = NULL,
    
    y =
      "Forecast Accuracy",
    
    caption =
      paste(
        "Temp accuracy: ±1°C.",
        "Wind speed accuracy: ±5 kt.",
        "Wind direction accuracy: ≤20°.",
        "Red diamonds = mean."
      )
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    
    plot.subtitle =
      element_text(
        size = 11
      ),
    
    axis.title =
      element_text(
        face = "bold"
      ),
    
    legend.position = "none"
  )


#######################################
## Cloud height and amount verification
cloud_audit <- noise_data[
  !is.na(cloud_base_error),
  .(
    n = .N,
    base_accuracy = mean(abs(cloud_base_error) <= 500),
    
    base_mae = mean(abs(cloud_base_error)),
    
    amount_accuracy = mean(abs(cloud_amount_error) <= 1),
    
    amount_mae = mean(abs(cloud_amount_error))
  ),
  by = ForecasterID
]

##### 
############################################
## Combine Forecast Accuracy Data
## Temperature, Wind Speed, Wind Direction,
## Cloud Base, Cloud Amount
############################################

# 1. Accuracy data for temperature, wind speed, wind direction
acc_long <- rbindlist(list(
  temp_audit[, .(ForecasterID, Variable = "Temperature", accuracy)],
  wspd_audit[, .(ForecasterID, Variable = "Wind Speed", accuracy)],
  wdir_audit[, .(ForecasterID, Variable = "Wind Direction", accuracy)],
  cloud_audit[, .(ForecasterID, Variable = "Cloud Base", accuracy = base_accuracy)],
  cloud_audit[, .(ForecasterID, Variable = "Cloud Amount", accuracy = amount_accuracy)]
),
fill = TRUE)

############################################
## 2. MAE data
############################################

mae_long <- melt(
  accuracy_metric,
  id.vars = "ForecasterID",
  measure.vars = c(
    "temp_mae",
    "wind_spd_mae",
    "wind_dir_mae",
    "cloud_base_mae",
    "cloud_amount_mae"
  ),
  variable.name = "mae_variable",
  value.name = "mae"
)

mae_long[, Variable := fcase(
  mae_variable == "temp_mae",
  "Temperature",
  mae_variable == "wind_spd_mae",
  "Wind Speed",
  mae_variable == "wind_dir_mae",
  "Wind Direction",
  mae_variable == "cloud_base_mae",
  "Cloud Base",
  mae_variable == "cloud_amount_mae",
  "Cloud Amount"
)]

mae_long <- mae_long[, .(ForecasterID, Variable, mae)]

############################################
## 3. Combined accuracy and MAE data
############################################
acc_combined <- merge(acc_long,
                      mae_long,
                      by = c("ForecasterID", "Variable"),
                      all = TRUE)

acc_combined[, Variable := factor(
  Variable,
  levels = c(
    "Temperature",
    "Wind Speed",
    "Wind Direction",
    "Cloud Base",
    "Cloud Amount"
  )
)]

############################################
## 4. Summary table for report
############################################
accuracy_summary <- acc_combined[, .(
  n_forecasters =
    uniqueN(ForecasterID),
  mean_accuracy =
    mean(accuracy, na.rm = TRUE),
  sd_accuracy =
    sd(accuracy, na.rm = TRUE),
  noise_index_pp =
    (2 / sqrt(pi)) *
    sd(accuracy, na.rm = TRUE) *
    100,
  mean_mae =
    mean(mae, na.rm = TRUE),
  sd_mae =
    sd(mae, na.rm = TRUE)
), by = Variable][order(Variable)]

accuracy_summary_report <- accuracy_summary[
  ,
  .(
    Variable,
    n_forecasters,
    mean_accuracy =
      round(mean_accuracy, 3),
    sd_accuracy =
      round(sd_accuracy, 3),
    noise_index_pp =
      round(noise_index_pp, 2),
    mean_mae =
      round(mean_mae, 2),
    sd_mae =
      round(sd_mae, 2)
  )
]

accuracy_summary_report

accuracy_summary_report_pct <- accuracy_summary[
  ,
  .(
    Variable,
    n_forecasters,
    mean_accuracy_pct =
      round(mean_accuracy * 100, 1),
    sd_accuracy_pp =
      round(sd_accuracy * 100, 1),
    noise_index_pp =
      round(noise_index_pp, 2),
    mean_mae =
      round(mean_mae, 2),
    sd_mae =
      round(sd_mae, 2)
  )
]

accuracy_summary_report_pct

ggplot(
  acc_combined,
  aes(
    x = Variable,
    y = accuracy,
    fill = Variable
  )
) +
  geom_violin(
    alpha = 0.4,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.size = 0.8
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "red",
    colour = "black"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 1)
  ) +
  labs(
    title = "Forecast Accuracy Noise",
    subtitle = "Forecaster accuracy distributions for temperature, wind, and cloud forecasts",
    x = NULL,
    y = "Forecast Accuracy",
    caption = paste(
      "Temperature: ±1°C;",
      "Wind speed: ±5 kt;",
      "Wind direction: ≤20°;",
      "Cloud base: ±500 ft;",
      "Cloud amount: ±1 category.",
      "Red diamonds = mean."
    )
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    plot.subtitle =
      element_text(
        size = 11
      ),
    axis.title =
      element_text(
        face = "bold"
      ),
    legend.position = "none"
  )

############################################
## Airport-adjusted forecaster accuracy
############################################
## Temperature
temp_acc_airport <- temp_data[, .(
  accuracy = mean(abs(temp_error) <= 1, na.rm = TRUE), n = .N), 
  by = .(ForecasterID, aerodrome)][, Variable := "Temperature"]

## Wind speed
wspd_acc_airport <- noise_data[, .(
  accuracy = mean(abs(wind_spd_error) <= 5, na.rm = TRUE), n = .N), 
  by = .(ForecasterID, aerodrome)][, Variable := "Wind Speed"]

## Wind direction
wdir_acc_airport <- noise_data[, .(
  accuracy = mean(abs(wind_dir_error) <= 20, na.rm = TRUE), n = .N), 
  by = .(ForecasterID, aerodrome)][, Variable := "Wind Direction"]

## Cloud base
cbase_acc_airport <- noise_data[!is.na(cloud_base_error), .(
  accuracy = mean(abs(cloud_base_error) <= 500, na.rm = TRUE), n = .N), 
  by = .(ForecasterID, aerodrome)][, Variable := "Cloud Base"]

## Cloud amount
camount_acc_airport <- noise_data[!is.na(cloud_amount_error), .(
  accuracy = mean(abs(cloud_amount_error) <= 1, na.rm = TRUE), n = .N), 
  by = .(ForecasterID, aerodrome)][, Variable := "Cloud Amount"]

## Combine
acc_airport <- rbindlist(
  list(
    temp_acc_airport,
    wspd_acc_airport,
    wdir_acc_airport,
    cbase_acc_airport,
    camount_acc_airport
  ),
  fill = TRUE
)

acc_airport[
  ,
  Variable := factor(
    Variable,
    levels = c(
      "Temperature",
      "Wind Speed",
      "Wind Direction",
      "Cloud Base",
      "Cloud Amount"
    )
  )
]

############################################
## Airport-adjusted accuracy anomalies
############################################

acc_airport[, airport_mean_accuracy :=
              weighted.mean(accuracy, w = n, na.rm = TRUE), by = .(aerodrome, Variable)]

acc_airport[, accuracy_anomaly :=
              accuracy - airport_mean_accuracy]

##
# accuracy_anomaly > 0 = better than peer average at that aerodrome for that variable
# accuracy_anomaly < 0 = below peer average at that aerodrome for that variable

############################################
## Forecaster-variable adjusted accuracy
############################################

forecaster_var_accuracy <- acc_airport[, .(
  adjusted_accuracy_anomaly =
    weighted.mean(accuracy_anomaly, w = n, na.rm = TRUE),
  
  total_n =
    sum(n, na.rm = TRUE)
), by = .(ForecasterID, Variable)]

############################################
## Standardised adjusted accuracy index
############################################

forecaster_var_accuracy[, adjusted_accuracy_z :=
                          as.numeric(scale(adjusted_accuracy_anomaly)), by = Variable]

##
# +1 = one standard deviation above peer-adjusted average
# 0 = average
# -1 = one standard deviation below peer-adjusted average

############################################
## Overall adjusted forecaster accuracy index
############################################

forecaster_accuracy_index <- forecaster_var_accuracy[, .(
  adjusted_accuracy_index =
    mean(adjusted_accuracy_z, na.rm = TRUE),
  
  n_variables =
    sum(!is.na(adjusted_accuracy_z)),
  
  mean_raw_anomaly_pp =
    mean(adjusted_accuracy_anomaly, na.rm = TRUE) * 100
), by = ForecasterID][order(-adjusted_accuracy_index)]

# Table
accuracy_differential_summary <- forecaster_accuracy_index[, .(
  n_forecasters = .N,
  mean_index =
    mean(adjusted_accuracy_index, na.rm = TRUE),
  sd_index =
    sd(adjusted_accuracy_index, na.rm = TRUE),
  min_index =
    min(adjusted_accuracy_index, na.rm = TRUE),
  max_index =
    max(adjusted_accuracy_index, na.rm = TRUE)
)]

accuracy_differential_summary

top_bottom_accuracy <- rbind(
  forecaster_accuracy_index[1:10],
  forecaster_accuracy_index[(.N - 9):.N]
)

top_bottom_accuracy

############################################
## Plot: airport-adjusted forecaster accuracy
############################################

plot_accuracy_index <- forecaster_accuracy_index[order(adjusted_accuracy_index)]

plot_accuracy_index[, Direction :=
                      fifelse(
                        adjusted_accuracy_index >= 0,
                        "Above peer-adjusted average",
                        "Below peer-adjusted average"
                      )]

# Create ordered factor
plot_accuracy_index[,ForecasterID_ordered := factor(ForecasterID, levels = ForecasterID)]
# Label every 5th forecaster only
label_every <- 5
y_labels <- levels(plot_accuracy_index$ForecasterID_ordered)
y_labels_sparse <- ifelse(seq_along(y_labels) %% label_every == 0, y_labels,"")
names(y_labels_sparse) <- y_labels

# Plot
ggplot(plot_accuracy_index,
       aes(
         x = reorder(ForecasterID, adjusted_accuracy_index),
         y = adjusted_accuracy_index,
         fill = Direction
       )) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = 0,
             linetype = 2,
             colour = "black") +
  scale_x_discrete(labels = y_labels_sparse) +
  scale_fill_manual(
    values = c(
      "Above peer-adjusted average" = "steelblue",
      "Below peer-adjusted average" = "firebrick"
    )
  ) +
  labs(
    title = "Airport-Adjusted Forecaster Accuracy Differentials",
    subtitle = "Composite accuracy anomaly after adjusting for aerodrome and forecast variable",
    x = "Forecaster ID",
    y = "Adjusted Accuracy Index",
    fill = NULL,
    caption = paste(
      "Positive index: above-average peer-adjusted accuracy.",
      "Negative index: below-average peer-adjusted accuracy."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title =
      element_text(face = "bold", size = 16),
    axis.title =
      element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

# heatmap
############################################
## Heatmap: adjusted accuracy by variable
############################################
ggplot(forecaster_var_accuracy,
       aes(
         x = Variable,
         y = reorder(ForecasterID, adjusted_accuracy_z),
         fill = adjusted_accuracy_z
       )) +
  geom_tile() +
  scale_fill_gradient2(
    low = "firebrick",
    mid = "white",
    high = "steelblue",
    midpoint = 0,
    name = "Adjusted\naccuracy"
  ) +
  labs(
    title = "Forecaster Accuracy Anomalies by Forecast Element",
    subtitle = "Accuracy relative to airport and variable peer averages",
    x = NULL,
    y = "Forecaster ID"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title =
          element_text(face = "bold", size = 16),
        axis.title =
          element_text(face = "bold"))


###################################################
# Part 7 — Full Noise Decomposition
# The challenge is that Kahneman's decomposition requires repeated judgments of the same case by different judges. We don't have that directly. We therefore need operational proxies.
# Proposed Aviation Noise Decomposition

### Accuracy Component ### 
accuracy_metric <- merge(temp_data[
    ,
    .(
      temp_mae = mean(abs(temp_error),na.rm = TRUE)),by = ForecasterID],
  noise_data[
    ,
    .(wind_spd_mae = mean(abs(wind_spd_error),na.rm = TRUE),
      wind_dir_mae = mean(abs(wind_dir_error),na.rm = TRUE),
      cloud_base_mae = mean(abs(cloud_base_error),na.rm = TRUE),
      cloud_amount_mae = mean(abs(cloud_amount_error),na.rm = TRUE)
    ),
    by = ForecasterID
  ],
  by = "ForecasterID",
  all = TRUE
)

accuracy_metric

### Hazard Component ### 
hazard_metric <- noise_data[
  ,
  .(
    fg_miss_rate =  mean(abs(fg_error),na.rm = TRUE),
    ts_miss_rate = mean(abs(ts_error),na.rm = TRUE),
    ra_miss_rate = mean(abs(ra_error),na.rm = TRUE),
    shra_miss_rate = mean(abs(shra_error),na.rm = TRUE)
  ),
  by = ForecasterID
]

hazard_metric

###### Communication Component #####
communication_metric <- merge(
  structure_entropy[
    ,
    .(ForecasterID,entropy)
  ],
  
  structure_noise[
    ,
    .(structure_sd = sd(pct)),
    by = ForecasterID
  ],
  by = "ForecasterID",
  all = TRUE
)

communication_metric

cloud_audit <- noise_data[
  !is.na(cloud_base_error),
  .(
    n = .N,
    cloud_base_accuracy = mean(abs(cloud_base_error) <= 500),
    cloud_amount_accuracy = mean(abs(cloud_amount_error) <= 1),
    cloud_base_mae = mean(abs(cloud_base_error)),
    cloud_amount_mae = mean(abs(cloud_amount_error))
  ),
  by = ForecasterID
]

wx_audit <- noise_data[
  ,
  .(
    ra_accuracy = mean(ra_error == 0,na.rm = TRUE),
    shra_accuracy = mean(shra_error == 0,na.rm = TRUE),
    dz_accuracy = mean(dz_error == 0,na.rm = TRUE),
    br_accuracy = mean(br_error == 0,na.rm = TRUE)
  ),
  by = ForecasterID
]

#################################
### THREE INDEX CONSTRUCTION
######

## PART 1 — Hazard Judgment Noise Index (HJNI)
# [forecaster-level hazard miss rates]
hazard_index <- noise_data[
  ,
  .(
    
    fg_noise =
      mean(abs(fg_error), na.rm=TRUE),
    
    ts_noise =
      mean(abs(ts_error), na.rm=TRUE),
    
    ra_noise =
      mean(abs(ra_error), na.rm=TRUE),
    
    shra_noise =
      mean(abs(shra_error), na.rm=TRUE),
    
    dz_noise =
      mean(abs(dz_error), na.rm=TRUE),
    
    br_noise =
      mean(abs(br_error), na.rm=TRUE)
    
  ),
  by=ForecasterID
]

# Standardise:
hazard_cols <- names(hazard_index)[-1]

for(col in hazard_cols){
  
  hazard_index[
    ,
    paste0(col,"_z") :=
      as.numeric(scale(get(col)))
  ]
  
}

# Create HJNI:
hz_z <- grep(
  "_z$",
  names(hazard_index),
  value=TRUE
)

hazard_index[
  ,
  HJNI :=
    rowMeans(
      .SD,
      na.rm=TRUE
    ),
  .SDcols=hz_z
]

####
## PART 2 — Forecast Accuracy Noise Index (FANI)
## [Skill variability]
fani <- merge(
  
  temp_data[
    ,
    .(
      temp_mae =
        mean(
          abs(temp_error),
          na.rm=TRUE
        )
    ),
    by=ForecasterID
  ],
  
  noise_data[
    ,
    .(
      
      wind_spd_mae =
        mean(
          abs(wind_spd_error),
          na.rm=TRUE
        ),
      
      wind_dir_mae =
        mean(
          abs(wind_dir_error),
          na.rm=TRUE
        ),
      
      cloud_base_mae =
        mean(
          abs(cloud_base_error),
          na.rm=TRUE
        ),
      
      cloud_amount_mae =
        mean(
          abs(cloud_amount_error),
          na.rm=TRUE
        )
      
    ),
    by=ForecasterID
  ],
  
  by="ForecasterID",
  all=TRUE
  
)

# Standardise:
fani_cols <- names(fani)[-1]

for(col in fani_cols){
  
  fani[
    ,
    paste0(col,"_z") :=
      as.numeric(scale(get(col)))
  ]
}

# Create FANI:
fani_z <- grep(
  "_z$",
  names(fani),
  value=TRUE
)

fani[
  ,
  FANI :=
    rowMeans(
      .SD,
      na.rm=TRUE
    ),
  .SDcols=fani_z
]

######
## PART 3 — Communication Noise Index (CNI)
# [Entropy + structure variability]
structure_noise_forecaster <- structure_noise[
  ,
  .(
    structure_sd =
      sd(pct)
  ),
  by=ForecasterID
]

cni <- merge(
  
  structure_entropy,
  
  structure_noise_forecaster,
  
  by="ForecasterID",
  all=TRUE
  
)

# Standardise. High entropy <> high noise...therefore invert entropy:
cni[
  ,
  entropy_inverted :=
    max(entropy,na.rm=TRUE) -
    entropy
]

# Standardise:
cni[
  ,
  entropy_z :=
    as.numeric(
      scale(
        entropy_inverted
      )
    )
]

cni[
  ,
  structure_z :=
    as.numeric(
      scale(
        structure_sd
      )
    )
]

# Create CNI:
cni[
  ,
  CNI :=
    rowMeans(
      .SD
    ),
  .SDcols=
    c(
      "entropy_z",
      "structure_z"
    )
]

######################
## PART 4 — Composite TAF Noise Index
ctni <- Reduce(
  
  function(x,y)
    merge(
      x,
      y,
      by="ForecasterID",
      all=TRUE
    ),
  
  list(
    
    hazard_index[
      ,
      .(
        ForecasterID,
        HJNI
      )
    ],
    
    fani[
      ,
      .(
        ForecasterID,
        FANI
      )
    ],
    
    cni[
      ,
      .(
        ForecasterID,
        CNI
      )
    ]
    
  )
  
)

ctni[
  ,
  CTNI :=
    (
      HJNI +
        FANI +
        CNI
    ) / 3
]

#############################
##### FIGURES

## Figure 9 — Three-Dimensional Noise Audit
plot_data <- melt(
  
  ctni,
  
  measure.vars=
    c(
      "HJNI",
      "FANI",
      "CNI"
    ),
  
  variable.name=
    "NoiseDimension",
  
  value.name=
    "NoiseIndex"
  
)

ggplot(plot_data,aes(NoiseDimension,NoiseIndex,fill=NoiseDimension)) +
  geom_violin(alpha=0.4) +
  geom_boxplot(width=0.12,fill="white") +
  stat_summary(
    fun=mean,
    geom="point",
    shape=23,
    size=3,
    fill="red"
  ) +
  labs(
    title = "Three-Dimensional Aviation Forecast Noise Audit",
    subtitle = "Forecaster hazard judgement, forecast accuracy and communication noise",
    x=NULL,
    y= "Standardised Noise Index",
    caption = "Higher values indicate greater forecaster-specific noise."
  ) +
  theme_minimal(base_size=13) +
  theme(
    legend.position="none",
    plot.title=
      element_text(
        face="bold",
        size=16
      )
  )

################################
#### TABLES
# Table 7. Integrated Aviation Forecast Noise Audit
table_noise_audit <- Reduce(
  
  function(x,y)
    merge(
      x,
      y,
      by = "ForecasterID",
      all = TRUE
    ),
  
  list(
    
    hazard_index[
      ,
      .(
        ForecasterID,
        HJNI
      )
    ],
    
    fani[
      ,
      .(
        ForecasterID,
        FANI
      )
    ],
    
    cni[
      ,
      .(
        ForecasterID,
        CNI
      )
    ]
    
  )
  
)

table_noise_audit[
  ,
  CTNI :=
    rowMeans(
      cbind(
        HJNI,
        FANI,
        CNI
      ),
      na.rm = TRUE
    )
]

table_noise_audit <- table_noise_audit[
  order(
    -CTNI
  )
]

table_noise_audit

table_noise_paper <- rbind(
  
  table_noise_audit[1:10],
  
  table_noise_audit[
    (.N-9):.N
  ]
  
)

table_noise_paper

noise_summary <- data.table(
  
  Dimension = c(
    "Hazard Judgement",
    "Forecast Accuracy",
    "Communication",
    "Composite TAF"
  ),
  
  Mean = c(
    
    mean(
      table_noise_audit$HJNI,
      na.rm = TRUE
    ),
    
    mean(
      table_noise_audit$FANI,
      na.rm = TRUE
    ),
    
    mean(
      table_noise_audit$CNI,
      na.rm = TRUE
    ),
    
    mean(
      table_noise_audit$CTNI,
      na.rm = TRUE
    )
    
  ),
  
  SD = c(
    
    sd(
      table_noise_audit$HJNI,
      na.rm = TRUE
    ),
    
    sd(
      table_noise_audit$FANI,
      na.rm = TRUE
    ),
    
    sd(
      table_noise_audit$CNI,
      na.rm = TRUE
    ),
    
    sd(
      table_noise_audit$CTNI,
      na.rm = TRUE
    )
    
  ),
  
  Min = c(
    
    min(
      table_noise_audit$HJNI,
      na.rm = TRUE
    ),
    
    min(
      table_noise_audit$FANI,
      na.rm = TRUE
    ),
    
    min(
      table_noise_audit$CNI,
      na.rm = TRUE
    ),
    
    min(
      table_noise_audit$CTNI,
      na.rm = TRUE
    )
    
  ),
  
  Max = c(
    
    max(
      table_noise_audit$HJNI,
      na.rm = TRUE
    ),
    
    max(
      table_noise_audit$FANI,
      na.rm = TRUE
    ),
    
    max(
      table_noise_audit$CNI,
      na.rm = TRUE
    ),
    
    max(
      table_noise_audit$CTNI,
      na.rm = TRUE
    )
    
  )
  
)

noise_summary

###
## Statistical Outputs for Results
noise_summary <- data.table(
  
  Index = c(
    "HJNI",
    "FANI",
    "CNI",
    "CTNI"
  ),
  
  Mean = c(
    mean(hazard_index$HJNI, na.rm=TRUE),
    mean(fani$FANI, na.rm=TRUE),
    mean(cni$CNI, na.rm=TRUE),
    mean(ctni$CTNI, na.rm=TRUE)
  ),
  
  SD = c(
    sd(hazard_index$HJNI, na.rm=TRUE),
    sd(fani$FANI, na.rm=TRUE),
    sd(cni$CNI, na.rm=TRUE),
    sd(ctni$CTNI, na.rm=TRUE)
  )
)

noise_summary

# Ensure CTNI is only calculated where all three components exist
ctni[
  ,
  CTNI := fifelse(
    !is.na(HJNI) & !is.na(FANI) & !is.na(CNI),
    (HJNI + FANI + CNI) / 3,
    NA_real_
  )
]

noise_summary <- data.table(
  Dimension = c(
    "Hazard Judgement",
    "Forecast Accuracy",
    "Communication",
    "Composite TAF"
  ),
  
  Mean = c(
    mean(ctni$HJNI, na.rm = TRUE),
    mean(ctni$FANI, na.rm = TRUE),
    mean(ctni$CNI, na.rm = TRUE),
    mean(ctni$CTNI, na.rm = TRUE)
  ),
  
  SD = c(
    sd(ctni$HJNI, na.rm = TRUE),
    sd(ctni$FANI, na.rm = TRUE),
    sd(ctni$CNI, na.rm = TRUE),
    sd(ctni$CTNI, na.rm = TRUE)
  ),
  
  N = c(
    sum(!is.na(ctni$HJNI)),
    sum(!is.na(ctni$FANI)),
    sum(!is.na(ctni$CNI)),
    sum(!is.na(ctni$CTNI))
  )
)

noise_summary


###################
## BIAS + NOISE ANALYSIS
# Step 1 — Estimate Bias
# Hazard Bias
hazard_bias <- noise_data[
  ,
  .(
    fg_bias = mean(fg_error, na.rm=TRUE),
    ts_bias = mean(ts_error, na.rm=TRUE),
    ra_bias = mean(ra_error, na.rm=TRUE),
    shra_bias = mean(shra_error, na.rm=TRUE)
  ),
  by = ForecasterID
]

hazard_bias

# Continuous Variable Bias 
# (already estimated)
continuous_bias <- noise_data[
  ,
  .(
    wind_spd_bias =
      mean(wind_spd_error, na.rm=TRUE),
    wind_dir_bias =
      mean(wind_dir_error, na.rm=TRUE),
    cloud_base_bias =
      mean(cloud_base_error, na.rm=TRUE),
    cloud_amount_bias =
      mean(cloud_amount_error, na.rm=TRUE)
  ),
  by = ForecasterID
]

continuous_bias

#########
# Construct Bias Index
bias_metrics <- merge(hazard_bias,continuous_bias,by = "ForecasterID",all = TRUE)

# Standardise
bias_cols <- names(bias_metrics)[-1]

for(col in bias_cols){
  bias_metrics[
    ,
    paste0(col,"_z") :=
      as.numeric(
        scale(
          abs(
            get(col)
          )
        )
      )
  ]
}

# Overall Bias Index
bias_z <- grep("_z$",names(bias_metrics),value = TRUE)
bias_metrics[,BiasIndex := rowMeans(.SD,na.rm = TRUE),.SDcols = bias_z]

# Results (Bias component of MSE)
summary(bias_metrics$BiasIndex)

#BiasIndex = rowMeans(abs(zscores),na.rm=TRUE)

###############################
## Now Define Total Noise
#########
ctni[,TotalNoise := rowMeans(cbind(HJNI,FANI,CNI),na.rm = TRUE)]

## Level Noise (different forecasters having different average tendencies)
LevelNoise <- mean(
  c(level_noise_summary$ts_noise,
    level_noise_summary$fg_noise,
    level_noise_summary$vis_noise,
    level_noise_summary$cld_noise)
)
LevelNoise

## Pattern Noise (forecasters responding differently to similar situations)
PatternNoise <- sd(cni$CNI,na.rm = TRUE)
PatternNoise

## Stable Pattern Noise (where the strongest weather-dependent disagreement exists)
StablePatternNoise <- mean(
  c(
    hazard_NI[hazard=="TS",NI],
    hazard_NI[hazard=="VIS",NI],
    hazard_NI[hazard=="CLD",NI]
  )
)

StablePatternNoise

## Occasion Noise (residual)
TotalNoise <- sd(ctni$TotalNoise,na.rm = TRUE)
TotalNoise

OccasionNoise <- sqrt(pmax(PatternNoise^2 - StablePatternNoise^2,0))
OccasionNoise

##############################
decomp <- data.table(
  Component = c("Bias","Level Noise","Stable Pattern Noise","Occasion Noise"),
  Value = c(mean(bias_metrics$BiasIndex,na.rm = TRUE),
    LevelNoise,
    StablePatternNoise,
    OccasionNoise)
)
decomp

decomp[,Variance := Value^2]
decomp

# library(data.table)
# library(ggplot2)
# 
# decomp <- data.table(Component = c("Bias",
#                                    "Level Noise",
#                                    "Stable Pattern Noise",
#                                    "Occasion Noise"),
#   Value = c(
#     0.02631781,
#     0.05607466,
#     0.07745267,
#     0.98510737
#   )
# )
# 
# decomp[,
#   Variance := Value^2
# ]

decomp[,end := cumsum(Variance)]
decomp[,start := shift(end,fill = 0)]
decomp[,xpos := .I]

decomp

ggplot(decomp) +
  geom_rect(aes(
    xmin = xpos - 0.4,
    xmax = xpos + 0.4,
    ymin = start,
    ymax = end,
    fill = Component
    )
    ) +
  geom_text(aes(
      x = xpos,
      y = end,
      label = round(Variance, 3)
    ),
    vjust = -0.5,
    size = 4
  ) +
  scale_x_continuous(
    breaks = seq_len(nrow(decomp)),
    labels = decomp$Component
  ) +
  labs(
    title = "Bias–Noise Decomposition of Aviation Forecast Variability",
    subtitle = "Variance contributions following the Kahneman framework",
    x = NULL,
    y = expression(
        "Variance Component ("*sigma^2*")"
      ),
    caption =
      paste(
        "Components shown on the variance scale.",
        "Occasion Noise dominates overall forecast variability."
      )) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 20,hjust = 1)
  ) + ylim(c(0,0.05))

decomp[,Contribution := 100 * Variance / sum(Variance)]

ggplot(
  decomp,
  aes(
    x = "",
    y = Contribution,
    fill = Component
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      label =
        paste0(
          round(
            Contribution,
            1
          ),
          "%"
        )
    ),
    position =
      position_stack(
        vjust = 0.5
      ),
    colour = "white",
    fontface = "bold"
  ) +
  labs(
    title =
      "Contribution of Bias and Noise Components",
    subtitle =
      "Relative contribution to total forecast variability",
    x = NULL,
    y =
      "Contribution (%)"
  ) +
  theme_minimal(base_size = 13)

##### REMOVE OCCASION NOISE
decomp_small <- data.table(
  Component = c(
    "Bias",
    "Level Noise",
    "Stable Pattern Noise"
  ),
  Value = c(
    mean(bias_metrics$BiasIndex, na.rm=TRUE),
    LevelNoise,
    StablePatternNoise
  )
)

ggplot(decomp_small,
  aes(reorder(Component, Value),
    Value,
    fill = Component)
  ) +
  geom_col(width = 0.7) +
  coord_flip() +
  geom_text(aes(label = round(Value, 3)),hjust = -0.2) +
  labs(title = "Measured Components of Forecast Variability",
    subtitle = "Bias and directly estimated noise components",
    x = NULL,
    y = "Magnitude",
    caption = "Occasion Noise was not directly measurable and is discussed qualitatively."
  ) +
  expand_limits(y = 0.1) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title =
      element_text(face="bold")
  )

#######################

# Step 2 — Estimate Total Noise
NoiseIndex = (HJNI + FANI + CNI) / 3

# Step 3 — Noise Decomposition
# Level Noise
LevelNoise <-
  mean(
    c(
      sd(hazard_airport$ts_anomaly),
      sd(hazard_airport$fg_anomaly),
      sd(hazard_airport$vis_anomaly),
      sd(hazard_airport$cld_anomaly)
    )
  )
LevelNoise # Forecaster "default settings"

# Pattern Noise
PatternNoise <- sd(cni$CNI,na.rm=TRUE)
PatternNoise

# Stable Pattern Noise (some situations consistently generate disagreement)
StablePatternNoise <-
  mean(c(hazard_NI[hazard=="TS",NI],
      hazard_NI[hazard=="VIS",NI],
      hazard_NI[hazard=="CLD",NI]
    )
  )

# Occasion Noise
# [Same forecaster, Same case, Different judgment]
OccasionNoise <- max(TotalNoise - LevelNoise - StablePatternNoise,0)

#### Figure
decomp <- data.table(
  Component = c(
    "Bias",
    "Level Noise",
    "Pattern Noise",
    "Stable Pattern Noise",
    "Occasion Noise"
  ),
  Value = c(
    Bias,
    LevelNoise,
    PatternNoise,
    StablePatternNoise,
    OccasionNoise)
)

ggplot(
  decomp,
  aes(
    x="Forecast Variability",
    y=Value,
    fill=Component
  )
) +
  geom_col()



###############
# LEVEL NOISE
# Some forecasters systematically forecast differently from others.
# Examples:
# 
# More thunderstorms
# More Low Visibility
# More Low Cloud
# More TEMPO groups
# We already use the airport-adjusted anomalies
# Between-forecaster variation after removing airport effects
level_noise <- mean(
  c(
    sd(hazard_airport$ts_anomaly),
    sd(hazard_airport$fg_anomaly),
    sd(hazard_airport$vis_anomaly),
    sd(hazard_airport$cld_anomaly)
  )
)

###########
# STABLE PATTERN NOISE
# Some situations consistently generate more disagreement than others
# e.g.
# Convective situations -> high disagreement
# Clear stable situations -> low disagreement
# Need to build a proxy -> use forecast complexity
# Create a complexity score:
noise_data[
  ,
  complexity :=
    as.numeric(ts_fcst) +
    as.numeric(fg_fcst) +
    as.numeric(vis_fcst) +
    as.numeric(cld1000_fcst)
]

stable_pattern_noise <- noise_data[
  ,
  .(
    ts_sd =
      sd(as.numeric(ts_fcst)),
    fg_sd =
      sd(as.numeric(fg_fcst)),
    vis_sd =
      sd(as.numeric(vis_fcst)),
    cld_sd =
      sd(as.numeric(cld1000_fcst))
  ),
  by = complexity
]

stable_pattern_noise_value <-
  mean(
    stable_pattern_noise$ts_sd,
    stable_pattern_noise$fg_sd,
    stable_pattern_noise$vis_sd,
    stable_pattern_noise$cld_sd
  )

stable_pattern_noise <-
  mean(
    c(
      level_noise_summary$ts_noise,
      level_noise_summary$vis_noise,
      level_noise_summary$cld_noise
    )
  )
stable_pattern_noise

##################
# PATTERN NOISE
# Forecasters respond differently to the same type of situation.
pattern_noise <- mean(structure_summary$sd_pct)

##################
# OCCASION NOISE
# The hardest
occasion_noise <- 1 -(level_noise + pattern_noise + stable_pattern_noise)

total_noise <-
  sd(
    c(
      hazard_airport$ts_anomaly,
      hazard_airport$fg_anomaly,
      hazard_airport$vis_anomaly,
      hazard_airport$cld_anomaly
    )
  )

occasion_noise <-
  max(
    total_noise -
      level_noise -
      pattern_noise -
      stable_pattern_noise,
    0
  )


noise_decomp <- data.table(
  
  Component = c(
    
    "Level Noise",
    
    "Pattern Noise",
    
    "Stable Pattern Noise",
    
    "Occasion Noise"
    
  ),
  
  Value = c(
    
    level_noise,
    
    pattern_noise,
    
    stable_pattern_noise,
    
    occasion_noise
    
  )
  
)

ggplot(
  noise_decomp,
  aes(
    reorder(Component,Value),
    Value,
    fill = Component
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    
    title =
      "Aviation Forecast Noise Decomposition",
    
    subtitle =
      "Estimated contributions of major noise components",
    
    x = NULL,
    
    y =
      "Relative Noise Magnitude",
    
    caption =
      paste(
        "Level noise reflects differences between forecasters;",
        "Pattern noise reflects differing responses to similar situations;",
        "Stable pattern noise reflects weather-regime-dependent disagreement;",
        "Occasion noise represents residual forecast variability."
      )
    
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "none",
    plot.title =
      element_text(face = "bold")
  )
