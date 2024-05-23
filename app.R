#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



# TODO : 
# - adjust plot size

# README ------------------------------------------------------------------
# USE R version 4.0.5


# Load packages -----------------------------------------------------------

library(shiny)
library(drc)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(DT)
library(gridExtra)
library(datamods)
library(plotly)

# Some functions ----------------------------------------------------------

# Function to generate an info icon with a tooltip.
generateInfoIcon <- function(title) {
    tags$i(
        class = "glyphicon glyphicon-info-sign",
        style = "color:#0072B2;",  # Blue color for the icon.
        title = title  # Tooltip title.
    )
}

# Function to extract model function names and concatenate with model name.
getMspec <- function(mod) {
    mod$fct$names %>% paste0(., collapse = ".") %>% 
        paste0(mod$fct$name, ":", .)
}

# function to compute model predictions with upper and lower CI
my_PM <- function(mod, tmpDf) {
    pm <-
        predict(mod,
                newdata = tmpDf,
                interval = "prediction",
                level = 0.95)
    tmpDf$p <- pm[, 1]
    tmpDf$pmin <- pm[, 2]
    tmpDf$pmax <- pm[, 3]
    tmpDf$pmin[tmpDf$pmin < 0] <-
        0  # replaces values below 0, otherwise confidence band stops at this point
    tmpDf$pmax[tmpDf$pmax > 100] <-
        100  # replaces values over 100, otherwise confidence band stops at this point
    return(tmpDf)
}



# Function to build a table with ECx values and confidence intervals.
buildECtable <- function(mod,
                         unit = "µg/L",
                         respLev = 50,
                         CIlevel = 0.95,
                         covFUN) {
    # Calculate ECx values and confidence intervals using the mod object.
    EC <- ED(
        mod,
        respLev = respLev,
        "delta",
        level = CIlevel,
        display = FALSE,
        vcov. = covFUN
    ) %>% as.data.frame(row.names = NULL)
    # Calculate and add confidence interval width to the table.
    EC[, paste0(CIlevel * 100, "%-CI")] <- EC$Upper - EC$Lower
    # Reorder and round columns in the table.
    EC <- EC[, c(1, 5, 2)] %>% round(., 2)
    # Add unit and model name to the table.
    EC$unit <- unit
    EC$model <- mod[["fct"]][["name"]]
    # Add lower and upper limits to the table, rounding to two decimal places.
    EC$lowerLimit <- ifelse(is.na(mod[["fct"]][["fixed"]][[2]]),
                            mod[["coefficients"]][["c:(Intercept)"]],
                            mod[["fct"]][["fixed"]][[2]]) %>% round(., 2)
    EC$upperLimit <- ifelse(is.na(mod[["fct"]][["fixed"]][[3]]),
                            mod[["coefficients"]][["d:(Intercept)"]],
                            mod[["fct"]][["fixed"]][[3]]) %>% round(., 2)
    # Modify ECx column names for readability.
    EC$ECx <- rownames(EC) %>% sub("^e:1:", "EC", .)
    # Add AIC values to the table, rounding to one decimal place.
    EC$AIC <- AIC(mod) %>% round(., 1)
    # Add a column for model specification, used for plotting.
    EC$Model <- getMspec(mod)
    # Reorder columns and rename some for clarity.
    EC <- EC[, c(8, 1, 2, 3, 4, 5, 6, 7, 9, 10)]
    names(EC)[names(EC) == 'lowerLimit'] <- 'lower limit'
    names(EC)[names(EC) == 'upperLimit'] <- 'upper limit'
    return(EC)  # Return the built EC table.
}


# Function to generate simulated response data based on concentration levels.
generate_response_data <- function(conc = rep(c(0, 5, 25, 125, 625, 3125), 5), mean_sd = 10, clip_min = 0, clip_max = 100) {
    # Calculate mean response based on a logistic function of concentration.
    mean_resp <- plogis(1.5 * (log(conc) - log(50))) * 100
    # Generate response data by adding normally distributed noise to the mean response.
    resp <- mean_resp + rnorm(length(mean_resp), mean = 0, sd = mean_sd)
    # Clip response values within specified minimum and maximum limits.
    resp <- pmin(pmax(resp, clip_min), clip_max)
    # Create a data frame with concentration and response columns.
    df <- data.frame(conc = conc, resp = resp)
    return(df)  # Return the generated response data.
}

example_data <- generate_response_data()

# UI ----------------------------------------------------------------------

ui <- fluidPage(titlePanel("ShinyECx"),
                sidebarLayout(
                    sidebarPanel(
                        fluidRow(
                            textInput("substance", label = NULL, placeholder = "Substance Name"),
                        ),
                        fluidRow(
                            actionButton("launch_modal", "Launch upload window"),
                            actionButton("example_data", "Generate Example Data")
                        ),
                        br(),
                        fluidRow(
                            column(width = 6,
                                   numericInput(
                                       "d",
                                       label = tags$span("Upper Limit",
                                                         generateInfoIcon(title = "Type in upper limit")),
                                       value = 100
                                   )),
                            column(width = 6,
                                   numericInput("c", label = "Lower limit", value = 0))
                        ),
                        
                        fluidRow(
                            column(width = 6,
                                   textInput("unit", 
                                             label = "Concentration unit", 
                                             value = "µg/L")
                            ),
                            column(width = 6,
                                   textInput("ec_values",
                                             label = "ECx values ",
                                             placeholder = 50,
                                             value = 50
                                   )
                            )),
                        
                        checkboxGroupInput(
                            "models",
                            label = "Choose Models",
                            choices = c(
                                "L.3",
                                "L.3 fixed",
                                "LL.3",
                                "LL.3 fixed",
                                "LL.4",
                                "LL.4 fixed",
                                "LL.5",
                                "LL.5 fixed",
                                "W1.3",
                                "W1.3 fixed",
                                "W2.3",
                                "W2.3 fixed",
                                "W1.4",
                                "W1.4 fixed",
                                "W2.4",
                                "W2.4 fixed"
                            ),
                            selected = "LL.3 fixed"
                        )
                        
                    ),
                    mainPanel(
                        width = 6,
                        tabsetPanel(
                            type = "tabs",
                            tabPanel("Table", 
                                     tableOutput("table")),
                            tabPanel(
                                "Plot",
                                plotOutput("doseResponsePlot", width = "600px", height = "400px"),
                                br(),
                                p(strong("Plot controls")),
                                column(
                                    width = 4,
                                    p("Plotwidth in pixel"),
                                    textInput("plotwidth",
                                              value = "600",
                                              label= NULL),
                                    p("Transparency"),
                                    sliderInput(
                                        "trans",
                                        NULL,
                                        min = 0,
                                        max = 1,
                                        value = 0.2
                                    ),
                                    p("F-factor"),
                                    sliderInput(
                                        "Ffactor",
                                        NULL,
                                        min = 1,
                                        max = 500,
                                        value = 5
                                    )
                                    
                                ),
                                column(
                                    width = 4,
                                    p("x-axis label"),
                                    textInput("xlabel", label = NULL, value = NULL),
                                    p("y-axis label"),
                                    textInput("ylabel", label = NULL, value = NULL),
                                    p("x-axis limit"),
                                    numericInput("xmax", label = NULL, value = 10)
                                ),
                                column(
                                    width = 4,
                                    p("log x-axis"),
                                    selectInput(
                                        "islog",
                                        label = NULL,
                                        choices = c("Yes", "No"),
                                        selected = "Yes"
                                    ),
                                    p("plot width"),
                                    numericInput("plot_width", label = NULL, value = 9),
                                    p("plot heigth"),
                                    numericInput("plot_height", label = NULL, value = 11),
                                    sliderInput(
                                        inputId = "opt.cex",
                                        label = "Point Size (cex)",
                                        min = 0,
                                        max = 2,
                                        step = 0.25,
                                        value = 1
                                    ),
                                    sliderInput(
                                        inputId = "opt.cexaxis",
                                        label = "Axis Text Size (cex.axis)",
                                        min = 0,
                                        max = 2,
                                        step = 0.25,
                                        value = 1
                                    ),
                                    downloadButton('export')
                                )
                            ),
                            tabPanel("EC values", DT::dataTableOutput("ECtable"))
                        )
                    )
                ))




# Server function ---------------------------------------------------------

server <- function(input, output, session) {
    
    observeEvent(input$example_data, {
        df <- generate_response_data()
        vals$data <- df
    })
    
    observeEvent(input$launch_modal, {
        req(input$from)
        import_modal(
            id = "myid",
            from = input$from,
            title = "Import data to be used in application"
        )
    })
    
    imported <- import_server("myid", return_class = "tbl_df")
    
    
    output$table <- renderTable({
        req(imported$data())
        imported$data()
    })
    
    
    
    vals <-
        reactiveValues(
            p1 = NULL,
            p2 = NULL,
            t1 = NULL,
            minc = NULL,
            minc2 = NULL,
            maxc = NULL,
            F0 = NULL,
            data = example_data
        )
    
    observeEvent(imported$data(), {
        vals$data <- imported$data()
    })
    
    observe({
        req(vals$data, input$Ffactor)
        df <- vals$data
        
        #vals$maxc <- max(df$conc)
        vals$maxc <- ifelse(is.numeric(input$maxc),
                            input$maxc,
                            max(df$conc))
        vals$minc <- df$conc[df$conc != 0] %>% sort() %>% .[1]
        vals$minc2 <- df$conc[df$conc != 0] %>% unique() %>% sort() %>% .[2]
        vals$F0 <- (vals$minc2 / vals$minc) * input$Ffactor
    })
    
    # Store the fitted models
    models <- reactiveValues()
    
    observeEvent(input$launch_modal, {
        # Clear the stored models when a new data file is uploaded
        models$values <- NULL
    })
    
    observeEvent(input$models, {
        req(vals$data)
        df <- vals$data 
        
        # Fit models if not already stored
        # Model parameters: b: Slope; c: lower limit; d: upper limit, e:ED50
        # 3: names = c("b", "d", "e")
        # 4: names = c("b", "c", "d", "e")
        # 5: names = c("b", "c", "d", "e", "f")
        for (model in input$models) {
            if (!model %in% names(models$values)) {
                fct <- switch(
                    model,
                    "L.3" = L.3(),
                    "L.3 fixed" = L.3(fixed = c(NA, input$d, NA)),
                    "LL.3" = LL.3(),
                    "LL.3 fixed" = LL.3(fixed = c(NA, input$d, NA)),
                    "LL.4" = LL.4(),
                    "LL.4 fixed" = LL.4(fixed = c(NA, input$c, input$d, NA)),
                    "LL.5" = LL.5(),
                    "LL.5 fixed" = LL.5(fixed = c(NA, input$c, input$d, NA, NA)),
                    "W1.3" = W1.3(),
                    "W1.3 fixed" = W1.3(fixed = c(NA, input$d, NA)),
                    "W2.3" = W2.3(),
                    "W2.3 fixed" = W2.3(fixed = c(NA, input$d, NA)),
                    "W1.4" = W1.4(),
                    "W1.4 fixed" = W1.4(fixed = c(NA, input$c, input$d, NA)),
                    "W2.4" = W2.4(),
                    "W2.4 fixed" = W2.4(fixed = c(NA, input$c, input$d, NA))
                )
                models$values[[model]] <-
                    drm(resp ~ conc, data = df, fct = fct)
            }
        }
    })
    
    
    output$table <- renderTable({
        vals$data
    })
    
    combinedECtable <- reactiveVal(NULL)
    combinedECtable_plot <- reactiveVal(NULL)
    
    observeEvent({
        input$models
        input$unit
        input$ec_values
    },
    {
        req(input$models, input$unit, input$ec_values, vals$data)
        
        ec_values <- as.numeric(strsplit(input$ec_values, ",\\s*")[[1]])
        
        
        
        # Call the buildECtable function to compute EC table for the selected model
        ECdata <- lapply(ec_values, function(respLev) {
            # Call the buildECtable function for each respLev
            lapply(input$models, function(model) {
                mod <- models$values[[model]]
                buildECtable(mod,
                             unit = input$unit,
                             respLev = respLev,
                             covFUN = vcov)
            })
        })
        
        # Combine the EC tables into a single data frame
        combinedECdata <- do.call(rbind, do.call(c, ECdata))
        combinedECdata_plot <- combinedECdata
        
        vals$t1 <- tableGrob(format(
            combinedECdata,
            core.just = "center",
            base_size = 10
        ),
        rows = NULL)
        
        # Update the reactive value with the combined EC table
        combinedECtable(as.data.frame(combinedECdata))
        combinedECtable_plot(as.data.frame(combinedECdata_plot))
    })
    
    output$ECtable <- DT::renderDataTable({
        datatable(
            combinedECtable(),
            rownames = FALSE,
            extensions = 'Buttons',
            options = list(dom = 'Bfrtip',
                           buttons = list('csv', 'excel', 'pdf'))
        )
    })
    
    
    
    
    
    #})
    
    output$doseResponsePlot <- renderPlot({
        req(input$models,
            input$d,
            input$c,
            input$unit,
            vals$data,
            vals)
        
        
        df <- vals$data #%>% select(input$concCol, input$respCol)
        colnames(df) <- c("conc", "resp")
        # Shift the 0 concentration values to avoid issues with log scale plotting
        df$conc0 <- df$conc
        
        
        df$conc0[df$conc0 == 0] <- vals$minc / vals$F0
        plotData <- lapply(input$models, function(model) {
            mod <- models$values[[model]]
            
            # Predicting values with confidence intervals for all generated dose levels
            
            # Generating new dose levels on the log scale as support for the line.

            newdata <-
                expand.grid(conc = exp(seq(log((vals$minc / vals$F0)
                ), log(vals$maxc * input$xmax), length = 200)))

            
            vis1 <-
                my_PM(mod, newdata)  # data points for model visualization 1
            
            # append model specifications
            
            vis1$Model <- getMspec(mod)
            vis1
            
            
        })
        
        # Combine all model data into a single data frame
        combinedData <- do.call(rbind, plotData)
        # some variable
        
        max_break <-
            ifelse(is.numeric(input$d),
                   input$d,
                   max(df$resp, na.rm = TRUE) * 1.1)
        xlabel <- ifelse(input$xlabel == "",
                         paste0("Concentration [", input$unit , "]"),
                         input$xlabel)
        ylabel <- ifelse(input$ylabel == "",
                         paste0("Embryos with effect [%]"),
                         input$xlabel)

        x_scale <- if (input$islog == "Yes") {
            scale_x_log10(breaks = 10^(-10:100),
                          minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9)))
        } else {
            scale_x_continuous(limits = c(0, input$xmax),
                               if (input$islog == "Yes") {
                                   annotation_logticks(sides = "b", scaled = TRUE)
                               } else {
                                   NULL
                               })
        }
        # x_scale <- if (input$islog == "Yes") {
        #   scale_x_log10(breaks = 10^(-10:10),
        #                 minor_breaks = rep(1:9, 21) * (10^rep(-10:10, each = 9)))
        # } else {
        #   scale_x_continuous(limits = c(0, input$xmax))
        # }
        
        # Plotting the combined data
        vals$p1 <- ggplot() +
            geom_point(data = df, aes(x = conc0, y = resp), color = "black") +
            scale_y_continuous(breaks=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))+
            geom_ribbon(
                data = combinedData,
                aes(
                    x = conc,
                    ymin = pmin,
                    ymax = pmax,
                    col = Model,
                    fill = Model
                ),
                alpha = input$trans,
                linetype = "blank"
            ) +
            geom_line(
                data = combinedData,
                aes(x = conc, y = p, col = Model),
                alpha = 0.8,
                linetype = "solid",
                size = 0.8
            ) +
            x_scale +
            ggtitle(input$substance) +
            xlab(xlabel) +
            ylab(ylabel) +
            guides(col = guide_legend(title = "Model")) +
            geom_segment(data = combinedECtable_plot(),
                         aes(x= Estimate, xend= Estimate, y=0, yend = input$d,
                             color = Model), linetype = "dashed"
            )+
            # geom_text_repel(
            #   data = combinedECtable_plot(),
            #   aes(
            #     x = Estimate,
            #     y = max_break,
            #     label = ECx,
            #     color = Model
            #   ),
            #   angle = 90,
            #   size = 4,
            #   point.size = NA,
        #   min.segment.length = Inf,
        #   ylim = c(max_break, Inf),
        #   box.padding = 0.1,
        #   max.overlaps = 1
        # )+
        geom_text(
            data = combinedECtable_plot(),
            aes(
                x = Estimate,
                y = max_break*1.1,
                label = ECx,
                color = Model
            ),
            angle = 90,
            size = 4,
            point.size = NA,
            min.segment.length = Inf,
            ylim = c(max_break*1*1, Inf),
            box.padding = 0.1,
            check_overlap = T
        )+
            theme_bw() +
            theme(
                legend.position = c(0.1, 0.7),
                plot.title = element_text(size = 20),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 10),
                plot.margin = margin(t = 50)
            )
        vals$p1
    })
    
    # Download handler for exporting the PDF file,
    # note that the file naming only works when running the app in a browser
    
    output$export <- downloadHandler(
        filename = function() {
            paste0("ShinyECx_", input$substance, ".pdf")
        },
        content = function(file) {
            g_tbl <- arrangeGrob(vals$t1, vals$p1,
                                 nrow=2, ncol=1, respect = TRUE, heights = rep(1, 2))
            ggsave(filename = file, plot = arrangeGrob(g_tbl), device = "pdf", width = input$plot_width, height = input$plot_height)
        }
    )
    
}

# Run app -----------------------------------------------------------------

shinyApp(ui, server)
