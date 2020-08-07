# Benjamin H Pepper
# BHPepper@gmail.com
# linkedin.com/in/BHPepper

library(shiny)
library(rlang)
library(pROC)
library(ggplot2)
library(fastDummies)

get_mod = function() {
  
  return(readRDS('integrated_mod.RDS'))
}

get_dat = function() {
  
  return(readRDS('integrated_dat.RDS'))
}

get_scores = function() {
  
  return(readRDS('integrated_scores.RDS'))
}

get_labels = function() {
  
  return(readRDS('integrated_labels.RDS'))
}

prob_f = function(mod, prototype, dat) {
  prototype_dummy = dummy_cols(rbind(prototype,dat), remove_selected_columns = T,
                    remove_first_dummy=T)[1,]
  preds = predict(mod$mod2, prototype_dummy, type='response')
  
  return(preds)
}

dat = get_dat()
mod = get_mod()
scores = get_scores()
labels = get_labels()

types = unlist(lapply(dat, class))
names = names(types)

prototype = dat[1,]

widgets = list()

for(i in 1:length(types)) {
  if (types[i] %in% c('numeric')) {
    if (names[i] == 'AGE_AT_DIAGNOSIS')
    {
      low = round(min(dat[,names[i]], na.rm=T),0)
      high = round(max(dat[,names[i]], na.rm=T),0)
      widgets[[i]] = sliderInput(inputId = names[i],
                                 h3(names[i]),
                                 min = low,
                                 max = high,
                                 value = low,
                                 step = 1)
    }
    else
    {
      low = round(min(dat[,names[i]], na.rm=T),2)
      high = round(max(dat[,names[i]], na.rm=T),2)
      widgets[[i]] = sliderInput(inputId = names[i],
                                 h3(names[i]),
                                 min = low,
                                 max = high,
                                 value = low,
                                 step = round((high - low)/50,2))
    }
  } else if (types[i] %in% c('character', 'factor')) {
    widgets[[i]] = selectInput(inputId = names[i],
                               h3(names[i]),
                               choices = unique(na.omit(dat[,names[i]])))
  } else if (types[i] %in% c('integer')) {
    low = min(dat[,names[i]], na.rm=T)
    high = max(dat[,names[i]], na.rm=T)
    widgets[[i]] = sliderInput(inputId = names[i],
                               h3(names[i]),
                               min = low,
                               max = high,
                               value = low,
                               step = 1)
  }
}

ui = fluidPage(
  titlePanel("Hormone Therapy Decision Support Tool"),
  
  sidebarLayout(
    exec('sidebarPanel', !!!widgets),
    mainPanel(
      h3(textOutput('results'),
        br(),
        plotOutput(outputId = "ROCPlot"))
    )  
  )
)

server = function(input, output) {
  output$ROCPlot = renderPlot({
    roc_obj = roc(labels, scores)
    ggroc(roc_obj) + theme_minimal() + labs(title = paste0('ROC Curve for the Model - AUC: ',
                                                           round(auc(roc_obj), 3)))
  })
  
  output$results = renderText({
    for(i in 1:length(types)) {
      prototype[,names[i]] = input[[names[i]]]
    }

    probability = prob_f(mod, prototype, dat)
    paste0('Predicted Probability of Hormone Therapy: ', round(probability, 3))
  })
}

shinyApp(ui = ui, server = server)
