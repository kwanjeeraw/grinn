mainPanel(width=12,
          fluidRow(column(12,
                          mainPanel(width=12,
                                    h3("setGrinnDb"),
                                    p("Set the graph database location for the currently working environment, see ",
                                      a(href='http://kwanjeeraw.github.io/grinn/setdb.html',target='_blank','here'),' for argument details.'
                                    )
                          )#end mainPanel
          )),
          wellPanel(
            h4("Current database location:"),
            fluidRow(
              column(12, verbatimTextOutput("currentdb"))
            ),
            h4("Input arguments:"), helpText("* required field"),
            fluidRow(
              column(6, textInput(inputId='dburl', label='url *', value=""))
            ),
            hr(),
            actionButton("submit","Submit")
          )
)