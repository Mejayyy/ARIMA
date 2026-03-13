library(shiny)
library(shinyjs)
library(DT)
cfs <- function(ts_data, lag_acf, lag_pacf) {

  par(mfrow   = c(1, 2),
      mar     = c(4, 4, 2, 1),  
      cex.lab = 1.3,           
      cex.axis= 1.1,           
      cex.main= 1.5)           
  

  acf(ts_data,
      main    = "ACF",
      lag.max = lag_acf)
  

  pacf(ts_data,
       main    = "PACF",
       lag.max = lag_pacf)
  

  par(mfrow = c(1, 1))
}


stationarity_check <- function(ts_data) {

  adf_result <- adf.test(ts_data, alternative = "stationary")
  

  kpss_level <- kpss.test(ts_data, null = "Level")
  

  kpss_trend <- kpss.test(ts_data, null = "Trend")
  

  print(adf_result)
  print(kpss_level)
  print(kpss_trend)
}

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Interactive Time Series Analysis"),

  actionButton("toggleSidebar", "Show/Hide Controls"),
  sidebarLayout(

    div(id = "sidebar",
        sidebarPanel(
          numericInput("frequency", "Frequency :", value = 5, min = 1, step = 1),
          numericInput("lag_max_acf", "Max lag for ACF:", value = 40, min = 1, step = 1),
          numericInput("lag_max_pacf", "Max lag for PACF:", value = 40, min = 1, step = 1),
          numericInput("diff",         "Diff (order):",     value = 1,   min = 1, step = 1),
          numericInput("max_order",    "Max order to test (for ARIMA(p,0,p) :", value = 4,   min = 1, step = 1),
          numericInput("chosen_order", "Chosen order for testing ( ARIMA(p,0,p-1), and so on):", 
                       value = 1, min = 1, step = 1),
          
          fluidRow(
            column(4, numericInput("p_final", "p:", value = 0, min = 0, step = 1)),
            column(4, numericInput("d_final", "d:", value = 1, min = 0, step = 1)),
            column(4, numericInput("q_final", "q:", value = 1, min = 0, step = 1))
          ),
          
          

        )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Time Series Original",
                 plotOutput("tsPlot")
        ),
        tabPanel(
          "ACF/PACF Original",
          plotOutput("acfPacfPlot")
        ),
        
        tabPanel(
          "Stationarity Original",

          verbatimTextOutput("stationarityResults")
        ),
        tabPanel(
          "Box–Cox for Original",

          verbatimTextOutput("boxcoxResults"),

          plotOutput("boxcoxPlot")
        ),
        tabPanel(
          "Diff ACF/PACF",
          plotOutput("diffSeriesPlot", height = "300px"),
          plotOutput("diffAcfPacfPlot", height = "500px")
        ),
        

        
        tabPanel(
          "Diff Stationarity",
          verbatimTextOutput("diffStationarityResults")
        ),
        
        tabPanel(
          "Results",
          h4("Model Selection on Differenced Series"),
          DTOutput("resultsTable")
        ),
        tabPanel(
          "Results (chosen order)",
          h4("Model Selection on Differenced Series"),
          DTOutput("resultsTable_2")
        ),
        tabPanel(
          "Diagnostics",
          h4("ARIMA Model Coefficients & σ²"),
          verbatimTextOutput("modelStats"),
          hr(),
          h4("tsdiag Plots"),
          plotOutput("tsdiagPlot", height = "400px"),
          hr(),
          h4("QQ Plot of Standardized Residuals"),
          plotOutput("qqPlot", height = "400px")
        ),
        tabPanel(
          "Forecast vs Actual",
          plotOutput("forecastPlot", height = "500px")
        )
        
        
        
        
        
        
        
      )
    )
  )

  
  )


server <- function(input, output, session) {

  observeEvent(input$toggleSidebar, {
    toggle("sidebar")
  })
  

  xau_df <- read.csv("XAU.csv", header = TRUE, stringsAsFactors = FALSE)
  close_all <- xau_df$Close.Last
  

  close_420 <- tail(close_all, 420)
  y_test    <- close_420[1:20]
  y_train   <- close_420[21:420]
  

  ts_train_reactive <- reactive({
    freq_val <- input$frequency
    ts(y_train, frequency = freq_val)
  })
  

  output$tsPlot <- renderPlot({
    ts_train <- ts_train_reactive()
    


    plot(
      ts_train,
      main = paste0("Training series (last 400 values of Close.Last)  |  freq = ", input$frequency),
      xlab = "Time",
      ylab = "Close.Last",
      type = "l",
      lwd = 2,
      col = "steelblue"
    )
    

    abline(h = 0, lty = 2, col = "gray")
  })
  
  output$acfPacfPlot <- renderPlot({
    ts_train <- ts_train_reactive()
    l1       <- input$lag_max_acf
    l2       <- input$lag_max_pacf
    

    cfs(ts_train, lag_acf = l1, lag_pacf = l2)
  })
  
  output$stationarityResults <- renderPrint({
    ts_train <- ts_train_reactive()
    stationarity_check(ts_train)
  })
  
  lambda_reactive <- reactive({
    ts_train <- ts_train_reactive()

    BoxCox.lambda(ts_train, lower = -3, upper = 3)
  })
  

  ts_train_bc_reactive <- reactive({
    ts_train <- ts_train_reactive()
    lam      <- lambda_reactive()
    BoxCox(ts_train, lam)
  })
  

  output$boxcoxResults <- renderPrint({
    lam    <- lambda_reactive()
    ts_bc  <- ts_train_bc_reactive()
    
    cat("Estimated lambda:", round(lam, 4), "\n\n")
    cat("Stationarity tests on Box–Cox–transformed series:\n")
    stationarity_check(ts_bc)
  })
  

  output$boxcoxPlot <- renderPlot({
    ts_bc <- ts_train_bc_reactive()
    lam   <- lambda_reactive()
    
    plot(
      ts_bc,
      main = paste0("Box–Cox (λ = ", round(lam, 4), ")"),
      xlab = "Time",
      ylab = expression(paste("Transformed ", x[t])),
      type = "l",
      lwd  = 2,
      col  = "darkgreen"
    )
    abline(h = 0, lty = 2, col = "gray")
  })
  
  
  ts_diff_reactive <- reactive({
    ts_orig <- ts_train_reactive()
    d_ord   <- input$diff
    if (d_ord < 1) {

      ts_orig
    } else {
      diff(ts_orig, differences = d_ord)
    }
  })
  
  

  output$diffSeriesPlot <- renderPlot({
    ts_diff <- ts_diff_reactive()
    d_ord   <- input$diff
    
    plot(
      ts_diff,
      main = if (d_ord >= 1) {
        bquote(paste(Delta^.(d_ord), " of ", ts[train]))
      } else {
        "No differencing applied (diff = 0)"
      },
      ylab = if (d_ord >= 1) expression(paste(Delta^.(d_ord), " x[t]")) else "x[t]",
      xlab = "Time",
      type = "l",
      lwd  = 2,
      col  = "blue"
    )
    abline(h = 0, lty = 2, col = "gray")
  })
  

  output$diffAcfPacfPlot <- renderPlot({
    ts_diff <- ts_diff_reactive()
    l1      <- input$lag_max_acf
    l2      <- input$lag_max_pacf
    

    cfs(ts_diff, lag_acf = l1, lag_pacf = l2)
  })
  

  output$diffStationarityResults <- renderPrint({
    ts_diff <- ts_diff_reactive()
    stationarity_check(ts_diff)
  })
  
  

  results_df_reactive <- reactive({
    ts_diff <- ts_diff_reactive()
    max_p   <- input$max_order
    

    orders <- lapply(1:max_p, function(i) c(i, 0, i))
    

    results <- data.frame(
      model = character(0),
      AIC   = numeric(0),
      BIC   = numeric(0),
      SS    = numeric(0),
      stringsAsFactors = FALSE
    )
    

    for (o in orders) {
      fit <- tryCatch(
        arima(ts_diff, order = o),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        results <- rbind(
          results,
          data.frame(
            model = paste0("(", paste(o, collapse = ","), ")"),
            AIC   = round(fit$aic, 3),
            BIC   = round(BIC(fit), 3),
            SS    = round(mean(fit$residuals^2), 4),
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    results
  })
  

  output$resultsTable <- renderDT({
    df <- results_df_reactive()
    datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = 5,
        lengthMenu = c(5, 10, 20),
        autoWidth  = TRUE
      )
    )
  })
  
  results_2_df_reactive <- reactive({
    ts_diff <- ts_diff_reactive()
    chosen_one   <- input$chosen_order
    
    
    orders <- list()
    for(i in 0:chosen_one) {
      orders[[  2*(i+1) -1  ]] <- c(i, 0, chosen_one)
      orders[[  2*(i+1)   ]] <- c(chosen_one, 0, i)
      
    }
    
  
    

    results <- data.frame(
      model = character(0),
      AIC   = numeric(0),
      BIC   = numeric(0),
      SS    = numeric(0),
      stringsAsFactors = FALSE
    )
    

    for (o in orders) {
      fit <- tryCatch(
        arima(ts_diff, order = o),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        results <- rbind(
          results,
          data.frame(
            model = paste0("(", paste(o, collapse = ","), ")"),
            AIC   = round(fit$aic, 3),
            BIC   = round(BIC(fit), 3),
            SS    = round(mean(fit$residuals^2), 4),
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    results
  })
  

  output$resultsTable_2 <- renderDT({
    df <- results_2_df_reactive()
    datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = 5,
        lengthMenu = c(5, 10, 20,30,40),
        autoWidth  = TRUE
      )
    )
  })
  
  model_reactive <- reactive({
    ts_train <- ts_train_reactive()
    p <- input$p_final
    d <- input$d_final
    q <- input$q_final
    

    if (is.null(p) || is.na(p) || p < 0 ||
        is.null(d) || is.na(d) || d < 0 ||
        is.null(q) || is.na(q) || q < 0) {
      return(NULL)
    }

    p <- as.integer(p)
    d <- as.integer(d)
    q <- as.integer(q)
    

    tryCatch(
      arima(ts_train, order = c(p, d, q)),
      error = function(e) NULL
    )
  })
  
  output$modelStats <- renderPrint({
    mod <- model_reactive()
    if (is.null(mod)) {
      cat("Model not fitted. Check (p,d,q) values.\n")
      return()
    }
    cat("Coefficients:\n")
    print(mod$coef)
    cat("\nSigma^2:\n")
    print(mod$sigma2)
  })
  
  output$tsdiagPlot <- renderPlot({
    mod <- model_reactive()
    if (is.null(mod)) {
      plot.new()
      title("No ARIMA model fitted; adjust (p,d,q).")
      return()
    }
    tsdiag(mod)
  })
  
  output$qqPlot <- renderPlot({
    mod <- model_reactive()
    if (is.null(mod)) {
      plot.new()
      title("No ARIMA model fitted; adjust (p,d,q).")
      return()
    }
    rs      <- mod$residuals
    stdres  <- rs / sqrt(mod$sigma2)
    qq      <- qqnorm(stdres, plot.it = FALSE)
    
    plot(
      qq$x, qq$y,
      main  = "QQ Plot of Standardized Residuals\nwith Regression Line",
      xlab  = "Theoretical Quantiles",
      ylab  = "Sample Quantiles",
      pch   = 20
    )

    fit <- lm(qq$y ~ qq$x)
    abline(fit, col = "magenta", lwd = 1)

    abline(0, 1, col = "blue", lty = 2)
    legend(
      "topleft",
      legend = c("Sample Quantiles", "Regression Line", "y = x"),
      col    = c("black", "magenta", "blue"),
      pch    = c(20,     NA,         NA),
      lty    = c(NA,     1,          2),
      lwd    = c(NA,     1,          1),
      bg     = "white"
    )
  })
  
  
  output$forecastPlot <- renderPlot({
    ts_train <- ts_train_reactive()
    

    model_0_1_1 <- model_reactive()

    auto_fit <- tryCatch(
      auto.arima(
        ts_train,
        max.p          = 3,
        max.q          = 3,
        max.P          = 3,
        max.Q          = 3,
        max.d          = 2,
        max.D          = 1,
        max.order      = 6,
        stepwise       = FALSE,
        approximation  = FALSE,
        trace          = FALSE
      ),
      error = function(e) NULL
    )
    

    if (is.null(model_0_1_1) || is.null(auto_fit)) {
      plot.new()
      title("Error: Could not fit one or both models.")
      return()
    }
    

    f011  <- predict(model_0_1_1, n.ahead = 20)
    fauto <- predict(auto_fit,      n.ahead = 20)
    
    pred011  <- as.numeric(f011$pred)
    predauto <- as.numeric(fauto$pred)
    

    actual <- c(y_train[1:20], y_test)
    

    plot(
      actual,
      type = "l",
      lwd  = 2,
      col  = "black",
      xlab = "Time (predictions start at 21)",
      ylab = "Value",
      main = "20‐Step Forecasts vs. Actuals"
    )
    lines(
      x   = 21:40,
      y   = pred011,
      col = "red",
      lwd = 2,
      lty = 1
    )
    lines(
      x   = 21:40,
      y   = predauto,
      col = "blue",
      lwd = 2,
      lty = 2
    )
    legend(
      "topleft",
      legend = c(
        "Actual (y_train[1:20] + y_test)",
        "ARIMA(0,1,1) Forecast",
        "auto.arima Forecast"
      ),
      col = c("black", "red", "blue"),
      lty = c(1,       1,      2),
      lwd = 2,
      bty = "n"
    )
  })
  
  
}


shinyApp(ui = ui, server = server)



