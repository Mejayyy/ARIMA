library(shiny)
library(shinyjs)


ui <- fluidPage(
  useShinyjs(),
  titlePanel("Interactive Time Series Analysis"),

  actionButton("toggleSidebar", "Show/Hide Controls"),
  verbatimTextOutput("detected_p"),
  sidebarLayout(

    div(id = "sidebar",
        sidebarPanel(
          numericInput("n_obs", "Number of observations:", value = 200, min = 10, step = 10),
          numericInput("seed", "Random seed:", value = 420, min = 1, step = 1),
          sliderInput("alpha", "Significance level (alpha):",
                      min = 0.01, max = 0.99, value = 0.05, step = 0.01),
          numericInput("lag_max_acf", "Max lag for ACF:", value = 40, min = 1, step = 1),
          numericInput("lag_max_pacf", "Max lag for PACF:", value = 40, min = 1, step = 1),
          numericInput("lag_max_emp_pacf", "Max lag for Empirical PACF:", value = 50, min = 1, step = 1),
          numericInput("p_order", "AR order (p):", value = 2, min = 0, step = 1),
          uiOutput("ar_coeffs")
        )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Time Series",
                 plotOutput("tsPlot")
        ),
        tabPanel("ACF",
                 plotOutput("acfPlot")
        ),
        tabPanel("PACF",
                 plotOutput("pacfPlot")
        ),
        tabPanel("My (Empirical) ACF",
                 plotOutput("emp_acf")
        ),
        tabPanel("R ACF minus My ACF",
                 plotOutput("acf_diff")
        ),
       tabPanel("My PACF",
                plotOutput("emp_pacf")
                          
        )
      )
    )
  )
)


server <- function(input, output, session) {

  observeEvent(input$toggleSidebar, {
    toggle("sidebar")
  })
  

  output$ar_coeffs <- renderUI({
    req(input$p_order)
    lapply(seq_len(input$p_order), function(i) {
      numericInput(
        inputId = paste0("ar", i),
        label   = paste0("AR coefficient phi[", i, "]: "),
        value   = switch(as.character(i),
                         "1" = 1,
                         "2" = -0.9,
                         "3" = 0.5,
                         "4" = -0.2,
                         0)
      )
    })
  })
  

  threshold <- reactive({
    qnorm(1 - input$alpha/2) / sqrt(input$n_obs)
  })
  

  tsData <- reactive({
    set.seed(input$seed)

    ar_coefs <- if (input$p_order > 0) {
      sapply(seq_len(input$p_order), function(i) input[[paste0("ar", i)]])
    } else {
      numeric(0)
    }
    x <- arima.sim(
      model = list(order = c(input$p_order, 0, 0), ar = ar_coefs),
      n     = input$n_obs
    )
    ts(x)
  })
  

  gamma_re <- reactive({
    x    <- as.numeric(tsData())
    n    <- input$n_obs
    xbar <- mean(x)
    dev  <- x - xbar
    cov_fun <- function(h) sum(dev[(h+1):n] * dev[1:(n-h)]) / n
    g <- sapply(0:(n-1), cov_fun)
    list(g0 = g[1], g  = g[-1])
  })
  

  output$tsPlot <- renderPlot({
    plot(tsData(),
         main = "Simulated Time Series",
         ylab = "Value",
         xlab = "Time")
  })
  

  output$acfPlot <- renderPlot({
    acf(tsData(),
        lag.max = input$lag_max_acf,
        main    = paste0("ACF (alpha = ", input$alpha, ")"),
        ci      = 1 - input$alpha)
  })
  

  output$pacfPlot <- renderPlot({
    pacf(tsData(),
         lag.max = input$lag_max_pacf,
         main    = paste0("PACF (alpha = ", input$alpha, ")"),
         ci      = 1 - input$alpha)
  })
  

  output$emp_acf <- renderPlot({
    gr <- gamma_re()
    rho <- c(gr$g0, gr$g) / gr$g0
    ylim <- range(c(rho, threshold(), -threshold()))
    plot(rho,
         type = "h",
         lwd  = 1,
         xlab = "Lag",
         ylab = "Empirical Autocorrelation",
         main = "My (Empirical) ACF",
         ylim = ylim)
    abline(h =  threshold(), col = "magenta", lty = 2)
    abline(h = -threshold(), col = "magenta", lty = 2)
    abline(h = 0,           col = "black")
  })
  

  output$acf_diff <- renderPlot({
    acf_res <- acf(tsData(),lag.max =input$n_obs,  plot = FALSE)
    acfR <- as.numeric(acf_res$acf)
    gr   <- gamma_re()
    rho  <- c(gr$g0, gr$g) / gr$g0
    diff <- acfR - rho
    plot(diff,
         type = "h",
         lwd  = 1,
         xlab = "Lag",
         ylab = "Diff",
         main = "R ACF - My ACF",
         ylim = range(c(diff, 0)))
    abline(h = 0, col = "black")
  })
  
  
  
  
  pacf_re <- reactive({
    gr <- gamma_re()
    gamma <- c(gr$g0, gr$g)
    n_obs <- input$n_obs
    pacf_emp <- numeric(n_obs - 1)
    phi_vec <- numeric(0)
    V <- gamma[1]
    for(k in 1:(n_obs - 1)){
      if(k == 1) {
        phi_nn <- gamma[2] / gamma[1]
        phi_vec <- phi_nn
        V <- (gamma[1]^2 - gamma[2]^2) / gamma[1]
      } else {
        phi_nn <- (gamma[k+1] - sum(phi_vec * rev(gamma[2:k]))) / V
        phi_vec <- c(phi_vec - phi_nn * rev(phi_vec), phi_nn)
        V <- V * (1 - phi_nn^2)
      }
      pacf_emp[k] <- phi_nn
    }
    pacf_emp
  })
  
  
  output$emp_pacf <- renderPlot({
    pacf_emp <- pacf_re()
    maxlag <- input$lag_max_emp_pacf
    threshold_val <- threshold()
    ylim <- c(min(-threshold_val, pacf_emp[1:maxlag]),
              max(threshold_val, pacf_emp[1:maxlag]))
    plot(pacf_emp[1:maxlag],
         type = "h",
         lwd  = 1,
         xlab = "Lag",
         ylab = "PACF",
         main = "My Empirical PACF",
         ylim = ylim)
    abline(h = 0, col = "black")
    abline(h =  threshold_val, col = "magenta", lty = 2)
    abline(h = -threshold_val, col = "magenta", lty = 2)
  })
  
  p_detect <- reactive({
    pacf_emp <- pacf_re()
    idxs <- which(abs(pacf_emp) >= threshold())
    p <- 1
    if (length(idxs) == 0) {
      p <- 0
    } else {
      for (i in seq(length(idxs), 1)) {
        if (i == 1) {
          p <- 1
          break
        }
        if (i == 2 && idxs[i] - 1 == idxs[i-1]) {
          p <- 2
          break
        }
        if (i > 2 && idxs[i] - 1 == idxs[i-1] && idxs[i-1] - 1 == idxs[i-2]) {
          p <- idxs[i]
          break
        }
      }
    }
    p
  })
  
  output$detected_p <- renderText({
    paste0("This AR proccess is of order ", p_detect())
  })
  
  
}


shinyApp(ui = ui, server = server)
