library(shiny)
library(deSolve)
library(dplyr)
library(shinyMatrix)
library(echarts4r)
library(DT)
library(tidyr)

#devtools::install_github('INWTlab/shiny-matrix') fixes #21 on their github

dummy_cohort_matrix <- matrix(c(24999, 1, 0, 2, 50000, 0, 0, 2, 25000, 0, 0, 2), nrow = 3, ncol = 4, byrow = TRUE,
                              dimnames = list(c("Young", "Adult", "Elderly"),
                                              c("Susceptible", "Infected", "Recovered", "Time to recover")))

dummy_contact_matrix <- matrix(rep(1, 12), nrow = 3, ncol = 3, byrow = TRUE,
                               dimnames = list(c("Young", "Adult", "Elderly"),
                                               c("Young", "Adult", "Elderly")))

dummy_infectivity_matrix <- matrix(c(rep(0.2, 5), rep(0.3, 7)), nrow = 3, ncol = 3, byrow = TRUE,
                                   dimnames = list(c("Young", "Adult", "Elderly"),
                                                   c("Young", "Adult", "Elderly")))

sir_model <- function(time, stocks, auxs, ...){
  
  with(as.list(c(stocks, auxs)),{ 
    
    states <- matrix(stocks, nrow = auxs$num_cohorts, ncol = auxs$num_stocks) # Convert stocks to a matrix
    
    # Auxiliaries
    # calculate aBeta by checking time
    aBeta = auxs$betas[[1]]
    for (i in 1:length(auxs$intervention_start)){
      if (time >= auxs$intervention_start[i] & time <= auxs$intervention_end[i]){
        aBeta = auxs$betas[[2]]
      }
    }
    
    aLambda <- aBeta %*% states[, 2] # infection rate based on beta and number currently infected
    
    # Flows
    fIR <- aLambda * states[, 1] # incidence rate based on the number susceptible and lambda
    fRR <- states[, 2] / auxs$delays # recovery rate based on number infected and recovery delay
    
    # Integral of Stocks
    d_sSusceptible_dt <- -fIR # remove the newly infected from stock
    d_sInfected_dt <- fIR - fRR # get the new number currently infected
    d_sRecovered_dt <- fRR # move those no longer infected to recovered
    
    # if(any(is.na(d_sSusceptible_dt))) browser() # enter function when NAs are found
    
    # Return a list of the important model values
    return(list(c(d_sSusceptible_dt,
                  d_sInfected_dt,
                  d_sRecovered_dt),
                IR = fIR,
                RR = fRR,
                Lambda = aLambda))
  })
}

ui <- fluidPage(
  titlePanel("[WIP] Modifiable SIR App"),
  sidebarPanel(
    width = 4,
    actionButton("browser", "browser"),
    #actionButton("btn_update_matrix", "update matrix"),
    tags$h4("Parameters"),
    sliderInput("ui_until", "Run model for x iterations", min = 1, max = 1095, step = 1, value = 365),
    sliderInput("ui_step", "Step size", min = 1, max = 7, step = 1, value = 0.125),
    sliderInput("ui_intervention_date", "Intervention between", min = 1, max = 2500, step = 1, value = c(200, 250)),
    sliderInput("ui_contact_intervention_effect", "Contact % Reduction", min = 0, max = 100, step = 1, value = 20, post = "%"),
    sliderInput("ui_infectivity_intervention_effect", "Infectivity % Reduction", min = 0, max = 100, step = 1, value = 20, post = "%"),
    tags$h5("Cohorts"),
    matrixInput(
      "Cohorts",
      value = dummy_cohort_matrix, class = "numeric",
      cols = list(names = TRUE, editableNames = FALSE, extend = FALSE), # TODO enable this to link to new stocks
      rows = list(names = TRUE, editableNames = TRUE, extend = TRUE)
    ),
    tags$h5("Contacts"),
    matrixInput("Contacts", value = dummy_contact_matrix, class = "numeric"),
    tags$h5("Infectivity"),
    matrixInput("Infectivity", value = dummy_infectivity_matrix, class = "numeric")
  ),
  mainPanel(
    width = 8,
    echarts4rOutput("plot"),
    DTOutput('table')
  )
)

server <- function(input, output, session) {
  
  observe({
    cohorts <- rownames(input$Cohorts) #c(rownames(input$Cohorts), "Test", "Test2")
    current <- input$Contacts
    new_cohorts <- setdiff(cohorts, rownames(current)) # find missing cohorts
    
    if(length(new_cohorts) > 0){
      old_length <- length(rownames(current))
      
      new_rows <- matrix(rep(0.1, length(rownames(current)) * length(new_cohorts)),
                         nrow = length(new_cohorts), dimnames = list(new_cohorts, rownames(current)))
      
      new_cols <- matrix(rep(0.1, length(cohorts) * length(new_cohorts)),
                         ncol = length(new_cohorts), dimnames = list(cohorts, new_cohorts))
      
      current <- rbind(current, new_rows)
      current <- cbind(current, new_cols)
      
      updateMatrixInput(session, "Contacts", current)
    }
  })
  
  observe({
    cohorts <- rownames(input$Cohorts) #c(rownames(input$Cohorts), "Test", "Test2")
    current <- input$Infectivity
    new_cohorts <- setdiff(cohorts, rownames(current)) # find missing cohorts
    
    if(length(new_cohorts) > 0){
      old_length <- length(rownames(current))
      
      new_rows <- matrix(rep(0.1, length(rownames(current)) * length(new_cohorts)),
                         nrow = length(new_cohorts), dimnames = list(new_cohorts, rownames(current)))
      
      new_cols <- matrix(rep(0.1, length(cohorts) * length(new_cohorts)),
                         ncol = length(new_cohorts), dimnames = list(cohorts, new_cohorts))
      
      current <- rbind(current, new_rows)
      current <- cbind(current, new_cols)
      
      updateMatrixInput(session, "Infectivity", current)
    }
  })
  
  model_results <- reactive({
    stock_columns <- c("Susceptible", "Infected", "Recovered")
    
    # Create a vector of time to run the model over
    time <- seq(0, input$ui_until, by = input$ui_step) # TODO figure out why get negative with short runs/large by
    
    # Define the number of population cohorts (e.g age groups) and stocks (e.g. S, I, and R)
    cohorts <- data.frame(input$Cohorts[!rownames(input$Cohorts) %in% c(""), ], stringsAsFactors = F, check.names = F)
    cohorts[] <- lapply(cohorts, as.integer)
    
    # Setup the initial stocks
    stocks <- c()
    for(i in stock_columns){
      stocks <- c(stocks, setNames(object = cohorts[, i], nm = paste(i, row.names(cohorts), sep = "_")))
    }
    
    population_totals <- rowSums(cohorts[, stock_columns])
    
    aEffectiveContactRate <- input$Contacts * input$Infectivity
    beta <- aEffectiveContactRate / population_totals
    
    aEffectiveContactRate2 <- (input$Contacts * ((100 - input$ui_contact_intervention_effect) / 100)) * (input$Infectivity * ((100 - input$ui_infectivity_intervention_effect) / 100))
    beta2 <- aEffectiveContactRate2 / population_totals
    
    betas <- list(b1 = beta, b2 = beta2)
    
    auxs <- list(delays = cohorts$`Time to recover`, num_cohorts = nrow(cohorts), num_stocks = length(stock_columns),
                 intervention_start = input$ui_intervention_date[1], intervention_end = input$ui_intervention_date[2], betas = betas)
    
    o <- as.data.frame(ode(y = stocks, times = time, func = sir_model, parms = auxs, method="euler"))
  })
  
  output$plot <- renderEcharts4r({
    model_results() %>%
      select(time, contains("nfect")) %>%
      pivot_longer(cols = -time, names_to = "cohort") %>%
      group_by(cohort) %>%
      e_chart(time) %>%
      e_line(value) %>%
      e_tooltip(trigger = "axis")
  })
  
  output$table <- renderDataTable({
    model_results() %>%
      datatable()
  })
  
  observeEvent(input$browser, {browser()})
}

shinyApp(ui, server)

