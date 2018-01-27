rm(list=ls())



# List of packages for session
.packages = c("tm", "slam", "dplyr",'readr',"ggplot2","LDAvis","Rmpfr","SnowballC","wordcloud",
              "data.table","topicmodels","knitr","stringi","RWeka","koRpus","shiny","visNetwork",
              "tidytext","tidyr","dplyr","lazyeval","stringr","lexRankr","DT","plotly","RColorBrewer")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

library(tm)
library(slam)
library(dplyr)
library(readr)
library(ggplot2)
library(LDAvis)
library(Rmpfr)
library(SnowballC)
library(wordcloud)
library(data.table)
library(topicmodels)
library(knitr)
library(Rmpfr)
library(stringi)
library(RWeka)
library(koRpus)
library(shiny)
library(data.table)
library(visNetwork)
library(tidytext)
library(tidyr)
library(dplyr)
library(lazyeval)
library(stringr)
library(lexRankr)
library(DT)
library(plotly)
library(RColorBrewer)


#################################### customizing stopwords ######################################
mystop<-c('completely', 'nothing', 'n/a', 'no', 'the','this','was','will','they','and','need',
          'two','got','though','around','get','want','asked','thats','only',
          'many','very','can','said','just','also','lot','one','when','back',
          'went','usual','much','use','ever','able','always','make', 'dont')


################################# DTM ################################################

# 
prepare_dtm <- function(data, enable_bigram, min, max, mystop) {
  result = list()
  if(enable_bigram) {
    data1 = str_replace_all(data,"[^[:graph:]]", " ") 
    corpus0 <- VCorpus(VectorSource(data1))
    corpus <- tm_map(tm_map(tm_map(tm_map(tm_map(corpus0, stripWhitespace), 
                                          content_transformer(tolower)), 
                                   removePunctuation),
                            removeNumbers),
                     removeWords, c(mystop, tm::stopwords("english"), tm::stopwords('SMART')))
    CorpusCopy <- corpus
    corpus <- tm_map(corpus, stemDocument)
    BigramTokenizer <- function(x) NGramTokenizer(x, Weka_control(min = min, max = max))
    dtm <- DocumentTermMatrix(corpus, control = list(tokenize=BigramTokenizer ,weighting=weightTf))
    result$corpus <- corpus
    result$corpus0 <- corpus0
    result$CorpusCopy <- CorpusCopy
    result$dtm <- dtm
    return(result)
  }
  else {
    Corpus<-Corpus(VectorSource(data))
    dtm <- DocumentTermMatrix(Corpus,
                              control = list(tolower = TRUE, 
                                             removeNumbers = TRUE, 
                                             removePunctuation = TRUE,
                                             stopwords=c(mystop, stopwords('english')),
                                             stemming = TRUE,wordLengths = c(2,Inf)))
  }
}
############################# Cut words based on TF-IDF #################################
filter_dtm <- function(dtm, cutoff) {
  term_tfidf <- tapply(dtm$v/row_sums(dtm)[dtm$i],
                       dtm$j, mean) * log2(nDocs(dtm)/col_sums(dtm > 0))
  cat(summary(term_tfidf))
  threshold <- summary(term_tfidf)[cutoff]
  dtm <- dtm[, term_tfidf >= threshold]
}

######################## Remove zero word documents #####################################
remove_zero <- function(dtm) {
  result = list()
  result$toRemove <- which(row_sums(dtm) == 0)
  result$dtm <- dtm[row_sums(dtm) > 0,]
  return(result)
}


####################### Check optimal number of topics ###############################
optimal_k <- function(dtm) {
  harmonicMean <- function(logLikelihoods, precision=2000L) {
    llMed <- median(logLikelihoods)
    as.double(llMed - log(mean(exp(-mpfr(logLikelihoods,
                                         prec = precision) + llMed))))
  }
  burnin = 2000
  iter = 2000
  keep = 50
  sequ <- seq(2, 20, 2)
  fitted_many <- lapply(sequ, function(k) LDA(dtm, k = k,
                                              method = "Gibbs",
                                              control = list(burnin = burnin, 
                                                             iter = iter, 
                                                             keep = keep) ))
  logLiks_many <- lapply(fitted_many, function(L) L@logLiks[-c(1:(burnin/keep))])
  hm_many <- sapply(logLiks_many, function(h) harmonicMean(h))
  # plot(sequ, hm_many,type="l")
  df <- data.frame(nbr_topics = sequ, log_likelihood = hm_many)
  plot_ly(df, x = ~nbr_topics, y = ~log_likelihood) %>%
    add_lines()
  # geom_line()
  # a <- qplot(sequ, hm_many, geom='line', xlab='Number of Topics', ylab='Log_Likelihood')
  # ggplotly(a)
}

############################## Run Topic Modeling ####################################

tm <- function(dtm, nbr_topic) {
  k = nbr_topic
  SEED = 19
  models <- list(
    CTM       = CTM(dtm, k = k, control = list(seed = SEED, var = list(tol = 10^-4), em = list(tol = 10^-3)))
  )
}

######################################### Visualization #################################
visualize <- function(corpus, dtm1, dtm2, modelname) {
  topicmodels_json_ldavis <- function(fitted, corpus, doc_term){
    
    # Find required quantities
    phi <- posterior(fitted)$terms %>% as.matrix
    theta <- posterior(fitted)$topics %>% as.matrix
    vocab <- colnames(phi)
    doc_length <- vector()
    for (i in 1:length(corpus)) {
      temp <- paste(corpus[[i]]$content, collapse = ' ')
      doc_length <- c(doc_length, stri_count(temp, regex = '\\S+'))
    }
    # temp_frequency <- inspect(doc_term)
    temp_frequency <- as.matrix(doc_term)
    freq_matrix <- data.frame(ST = colnames(temp_frequency),
                              Freq = colSums(temp_frequency))
    rm(temp_frequency)
    
    # Convert to json
    json_lda <- LDAvis::createJSON(phi = phi, theta = theta,
                                   vocab = vocab,
                                   doc.length = doc_length,
                                   term.frequency = freq_matrix$Freq)
    
    return(json_lda)
  }
  Corpus_mod<-corpus[-remove_zero(dtm1)$toRemove]
  json_lda<-topicmodels_json_ldavis(modelname, Corpus_mod,dtm2)
  serVis(json_lda)
  # serVis(json_lda, out.dir = './', open.browser = FALSE)
}


############################# wordcloud ##########################################
prepare_wordcloud <- function(modelname, dtm, topic_nbr, filename) {
  topic <- topics(modelname, 1)
  ids <- as.numeric(dtm[as.numeric(names(topic[topic == topic_nbr])),]$dimnames$Docs)
  dtm1 <- dtm[ids, ]
  freq = sort(colSums(as.matrix(dtm1)),decreasing = TRUE)
  freq.df = data.frame(word=names(freq), freq=freq)
  pal=brewer.pal(8,"Blues")
  pal=pal[-(1:3)]
  # png(filename=paste(filename, "_topic_", topic_nbr, ".png", sep = ''))
  wordcloud(freq.df$word,freq.df$freq,max.words=100, random.order = F, colors=pal)
  # dev.off()
}

########################## Get original docs based on words ####################################
get_original_doc <- function(raw_data, dtm, phrase) {
  ids <- as.numeric(dtm[(as.matrix(dtm)[, phrase] != 0),]$dimnames$Docs)
  docs <- raw_data[ids,]
}


# 
############################### PUT low Christmas #######################################
run_tm <- function(raw_data, min = 2, max = 3, cutoff = 3, sparsity = 0.99, nbr_topics = 6, stopwords) {
  result = list()
  cat('prepare dtm')
  prepared = prepare_dtm(raw_data, enable_bigram = TRUE, min, max, stopwords)
  result$corpus <- prepared$corpus
  dtm <- prepared$dtm
  cat('cleaning dtm')
  dtm1 <- filter_dtm(dtm, cutoff)
  dtm2 <- removeSparseTerms(dtm1, sparsity)
  dtm3 <- remove_zero(dtm2)$dtm
  result$dtm1 <- dtm2
  result$dtm2 <- dtm3
  cat('modeling')
  result <- tm(dtm3, nbr_topics)
  # visualize(PUT_corpus, PUT_low_Chris_dtm2, PUT_low_Chris_dtm3, models$Gibbs)
  return(result)
}


############################# Prepare topics #####################################
prepare <- function(dtm, model, topic) {
  topicsProb <- topics(model, 1)
  dtm_mat<-as.matrix(dtm)
  topic_mat<-as.matrix(topicsProb)
  colnames(topic_mat)<-c("Topic")
  final_mat<-cbind(dtm_mat,topic_mat)
  print(paste0("topic = ", topic))
  final_mat_mod<-final_mat[final_mat[,"Topic"]==topic,]
  print("topic retrived...")
  Topic_data<-as.data.frame(final_mat_mod)
  Topic_data<-Topic_data[row_sums(Topic_data)>0,colSums(Topic_data) > 0]
  Topic_data[,"Topic"]<-NULL
  Topic_data$rownum <- as.numeric(rownames(Topic_data))
  
  Topic_data_tr<-melt(Topic_data, id = c("rownum"))
  Topic_data_tr<-filter(Topic_data_tr,Topic_data_tr$value>0)
  
  Topic_data_tr$variable<-as.character(Topic_data_tr$variable)
  print("get probability...")
  probability = posterior(model)$terms
  col_sum = col_sums(probability)
  prob = apply(probability, 1, function(x) x / col_sum)
  prob_topic = data.frame(prob[,topic])
  prob_topic$names = rownames(prob_topic)
  colnames(prob_topic) <- c('probability', 'words')
  rbPal <- colorRampPalette(c('red','pink'))
  prob_topic$color <- rbPal(10)[as.numeric(cut(prob_topic$probability,breaks = 10))]
  result = list()
  result$topic = Topic_data_tr
  result$prob = prob_topic
  return(result)
}

prepare_general <- function(dtm, model, topic) {
  topicsProb <- topics(model, 1)
  dtm_mat<-as.matrix(dtm)
  topic_mat<-as.matrix(topicsProb)
  colnames(topic_mat)<-c("Topic")
  final_mat<-cbind(dtm_mat,topic_mat)
  print(paste0("topic = ", topic))
  final_mat_mod<-final_mat[final_mat[,"Topic"]==topic,]
  print("topic retrived...")
  Topic_data<-as.data.frame(final_mat_mod)
  Topic_data<-Topic_data[row_sums(Topic_data)>0,colSums(Topic_data) > 0]
  Topic_data[,"Topic"]<-NULL
  Topic_data$rownum <- as.numeric(rownames(Topic_data))
  
  Topic_data_tr<-melt(Topic_data, id = c("rownum"))
  Topic_data_tr<-filter(Topic_data_tr,Topic_data_tr$value>0)
  
  Topic_data_tr$variable<-as.character(Topic_data_tr$variable)
  print("get probability...")
  probability = posterior(model)$terms
  col_sum = col_sums(probability)
  prob = apply(probability, 1, function(x) x / col_sum)
  prob_topic = data.frame(prob[,topic])
  prob_topic$names = rownames(prob_topic)
  colnames(prob_topic) <- c('probability', 'words')
  rbPal <- colorRampPalette(c('red','pink'))
  prob_topic$color <- rbPal(10)[as.numeric(cut(prob_topic$probability,breaks = 10))]
  result = list()
  result$topic = Topic_data_tr
  result$prob = prob_topic
  return(result)
}

############################# Create nodes and edges ######################################


create_node <- function(Topic_data_tr, prob, threshold) {
  nodes<-Topic_data_tr%>%
    select(-rownum)%>%
    group_by(variable)%>%
    summarise(num_word=sum(value))%>%
    arrange(desc(num_word))
  max_n <- max(nodes$num_word)
  if (max_n > threshold) {
    nodes <- nodes %>%
      filter(num_word >= threshold)
  }
  # filter(num_word >= threshold)
  colnames(nodes)<-c("id","size")
  nodes$size = nodes$size / 10.0
  
  nodes1 <- merge(prob, nodes, by.x = 'words', by.y = 'id')
  nodes1$label <- nodes1$words
  nodes1$id <- nodes1$words
  nodes1$words <- NULL
  return(nodes1)
}

create_edge <- function(Topic_data_tr, nodes) {
  temp1 <- Topic_data_tr %>% 
    select(-value) %>%
    filter(variable %in% nodes$id)
  
  temp2 <- Topic_data_tr %>% 
    select(-value) %>%
    filter(variable %in% nodes$id)
  
  join <- inner_join(temp1,temp2, by = "rownum")%>%
    filter(variable.x!=variable.y)
  
  edges<-join%>%
    group_by(from = variable.x,to = variable.y)%>%
    summarise(width=n() / 30.0)%>%
    arrange(desc(width))
  ids <- c(1:((dim(edges)[1])/2)) * 2
  edges <- edges[ids,]
  edges <- edges[1:50, ]
}

compute_data <- function(updateProgress = NULL) {
  # Create 0-row data frame which will be used to store data
  dat <- data.frame(x = numeric(0), y = numeric(0))
  
  for (i in 1:10) {
    Sys.sleep(0.25)
    
    # Compute new row of data
    new_row <- data.frame(x = rnorm(1), y = rnorm(1))
    
    # If we were passed a progress update function, call it
    if (is.function(updateProgress)) {
      text <- paste0("x:", round(new_row$x, 2), " y:", round(new_row$y, 2))
      updateProgress(detail = text)
    }
    
    # Add the new row of data
    dat <- rbind(dat, new_row)
  }
  
  dat
}


################################ prepare topic, node and edges ################################
topic_mat <- function(docs, dtm, modelname, topic_nbr) {
  topic <- topics(modelname, 1)
  ids <- as.numeric(dtm[as.numeric(names(topic[topic == topic_nbr])),]$dimnames$Docs)
  docs <- docs[ids, ]
  topic_dtm <- as.data.frame(as.matrix(prepare_dtm(docs$Openend, TRUE, 1, 1)$dtm))
  topic_dtm[row_sums(topic_dtm) <= 0,colSums(topic_dtm) <= 0]
  topic_dtm$rownum <- as.numeric(rownames(topic_dtm))
  topic_melt<-melt(topic_dtm, id = c("rownum"))
  topic_melt <- filter(topic_melt, topic_melt$value > 0)
  topic_melt$variable<-as.character(topic_melt$variable)
  return(topic_melt)
}



prepare_node <- function(topic_melt, threshold) {
  nodes <- topic_melt %>%
    select(-rownum) %>%
    group_by(variable) %>%
    summarise(num_word = sum(value)) %>%
    arrange(desc(num_word))%>%
    filter(num_word > threshold)
  colnames(nodes)<-c("id","size")
  nodes$size = nodes$size / 50.0
  return(nodes)
}

shape_to_tidy <- function(data, text, nbr_grams) {
  if (nbr_grams == 1) {
    tidy <- data %>%
      mutate(id = (1:dim(data)[1])) %>%
      mutate_(openend = interp(quote(as.character(var)), var = as.name(text))) %>%
      # select(id, openend) %>%
      unnest_tokens(word, openend, drop = TRUE) %>%
      filter(str_detect(word, "^[a-z]"),
             str_detect(word, "[a-z]$"),
             !word %in% stop_words$word) %>%
      select_(interp(quote(-(var)), var = as.name(text)))
  }
  if (nbr_grams == 2) {
    tidy <- data %>%
      mutate(id = (1:dim(data)[1])) %>%
      mutate_(openend = interp(quote(as.character(var)), var = as.name(text))) %>%
      # select(id, openend) %>%
      unnest_tokens(word, openend, token = "ngrams", n = 2) %>%
      separate(word, c("word1", "word2"), sep = " ") %>%
      filter(!word1 %in% stop_words$word) %>%
      filter(!word2 %in% stop_words$word) %>%
      unite(word, word1, word2, sep = " ") %>%
      select_(interp(quote(-(var)), var = as.name(text)))
  }
  return(tidy)
}

shape_to_tfidf <- function(tidydata, group_col, threshold) {
  tf_idf <- tidydata %>%
    mutate_(group = interp(quote(as.factor(var)), var = as.name(group_col))) %>%
    select(group, word) %>%
    group_by(group) %>%
    count(word, sort = TRUE) %>%
    arrange(group) %>%
    ungroup() %>%
    bind_tf_idf(word, group, n)
  max_n <- max(tf_idf$n)
  print("checking threshold status")
  print(max_n)
  print(threshold)
  if (threshold > max_n) {
    print("threshold larger than max frequency")
    threshold = n}
  tf_idf1 <- tf_idf %>%
    arrange(desc(tf_idf)) %>%
    group_by(group) %>%
    mutate(rank = row_number(tf_idf)) %>%
    filter(n >= threshold) %>%
    ungroup()
  return(tf_idf1)
}


nodes <- data.frame(id = 1:10,
                    
                    # add labels on nodes
                    label = paste("Node", 1:10),
                    
                    # add groups on nodes 
                    group = c("GrA", "GrB"),
                    
                    # size adding value
                    value = 1:10,          
                    
                    # control shape of nodes
                    shape = c("square", "triangle", "box", "circle", "dot", "star",
                              "ellipse", "database", "text", "diamond"),
                    
                    # tooltip (html or character), when the mouse is above
                    title = paste0("<p><b>", 1:10,"</b><br>Node !</p>"),
                    
                    # color
                    color = c("darkred", "grey", "orange", "darkblue", "purple"),
                    
                    # shadow
                    shadow = c(FALSE, TRUE, FALSE, TRUE, TRUE))             


edges <- data.frame(from = c(1,2,5,7,8,10), to = c(9,3,1,6,4,7))














server <- function(input, output, session) {
  options(shiny.maxRequestSize=60*1024^2)
  
  ################## Load data #################
  data <- reactive({
    cat("uploading files...")
    cat("/")
    if (is.null(input$file1)) {return()}
    file1 <- input$file1
    # read.table(file1$datapath, sep=input$sep, fileEncoding="UTF-8", fill = TRUE)
    data.frame(read_delim(file1$datapath, delim = input$sep1,
               escape_double = FALSE,
               # col_types = cols(Pickup_Date = col_date(format = "%Y-%m-%d"),                                                                                                                                         Survey_Completed_Date = col_date(format = "%Y-%m-%d")),
               locale = locale(encoding = "UTF-8")))
    # read_delim(file1$datapath, delim = input$sep1)
  })
  
  stopwords <- reactive({
    print('uploading stopwords...')
    mystop<-c('completely', 'nothing', 'n/a', 'no', 'the','this','was','will','they','and','need',
              'two','got','though','around','get','want','asked','thats','only',
              'many','very','can','said','just','also','lot','one','when','back',
              'went','usual','much','use','ever','able','always','make', 'dont')
    if (is.null(input$file2)) {return(mystop)}
    file2 <- input$file2
    read.delim(file2$datapath, sep = ',', header = FALSE)$V1
  })
  
  
  ############################ Topic Modeling #############################
  # run topic modeling based on user determined topic number
  
  dtm <- reactive({
    if (is.null(input$comment)) {return()}
    print(input$gram_range[1])
    print(input$gram_range[2])
    text <- data()[, c(input$comment)]
    print('preparing dtm...')
    prepared <- prepare_dtm(text, enable_bigram = TRUE, input$gram_range[1], input$gram_range[2], stopwords())
    dtm <- prepared$dtm
    corpus <- prepared$corpus0
    print('cleaning dtm...')
    dtm1 <- filter_dtm(dtm, cutoff = 2)
    dtm2 <- removeSparseTerms(dtm1, 0.99)
    remove <- remove_zero(dtm2)
    dtm3 <- remove$dtm
    toRemove <- remove$toRemove
    corpus_mod <- corpus[-toRemove]
    resp<-data.frame(text=unlist(sapply(corpus_mod, `[`, "content")), stringsAsFactors=F)
    dtms <- list()
    dtms$remove <- toRemove
    dtms$dtm1 <- dtm2
    dtms$dtm2 <- dtm3
    dtms$resp <- resp
    return(dtms)
  })
  
  result <- eventReactive(input$goButton, {
    print("running topic modeling...")
    
    progress <- shiny::Progress$new()
    progress$set(message = "Computing data", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Create a closure to update progress.
    # Each time this is called:
    # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
    #   distance. If non-NULL, it will set the progress to that value.
    # - It also accepts optional detail text.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    # Compute the new data, and pass in the updateProgress function so
    # that it can update the progress indicator.
    compute_data(updateProgress)
    
    if (is.null(input$nbr_topics)) {return()}
    n <- input$nbr_topics
    print(paste0("number of topics = ", n))
    result <- tm(dtm()$dtm2, nbr_topic = n)
    return(result)
  })
  
  data_for_download <- reactive({
    topicsProb <- topics(result()$CTM, 1)
    data1 <- data()[-dtm()$remove,]
    data1$topic <- topicsProb
    return(data1)
  })
  
  
  
  # Generate topics
  topic <- reactive({
    cat("preparing topics...")
    
    progress <- shiny::Progress$new()
    progress$set(message = "Generating topics...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Create a closure to update progress.
    # Each time this is called:
    # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
    #   distance. If non-NULL, it will set the progress to that value.
    # - It also accepts optional detail text.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    # Compute the new data, and pass in the updateProgress function so
    # that it can update the progress indicator.
    compute_data(updateProgress)
    
    if (is.null(input$topics)) {return()}
    n <- input$nbr_topics
    k <- as.numeric(strsplit(as.character(input$topics), "_")[[1]][2])
    topic <- prepare(dtm()$dtm2, result()$CTM, k)
    return(topic)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$file1$name, '.csv', sep = '')
    },
    content = function(file) {
      write.csv(data_for_download(), file)
    }
  )
  
  output$pie_chart <- renderPlotly({
    topicsProb <- topics(result()$CTM, 1)
    data <- data.frame(table(topicsProb))
    colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)', 'rgb(114,147,203)')
    
    p <- plot_ly(data, labels = ~paste("topic ", topicsProb, sep = ' '), values = ~Freq, type = 'pie',
                 textposition = 'inside',
                 textinfo = 'label+percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'text',
                 text = ~paste('frequency = ', Freq),
                 marker = list(colors = colors,
                               line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = FALSE) %>%
      layout(title = 'Topic Distribution Chart',
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    p
    
  })
  
  output$top_sentence <- renderTable({
    print("Generating the sentence...")
    
    progress <- shiny::Progress$new()
    progress$set(message = "Generating top sentences...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Create a closure to update progress.
    # Each time this is called:
    # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
    #   distance. If non-NULL, it will set the progress to that value.
    # - It also accepts optional detail text.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    # Compute the new data, and pass in the updateProgress function so
    # that it can update the progress indicator.
    compute_data(updateProgress)
    
    
    if (is.null(input$topics1)) {return()}
    if (is.null(input$s_threshold)) {return()}
    if (is.null(input$s_nbr)) {return()}
    topicsProb <- topics(result()$CTM, 1)
    
    
    resp<-cbind(dtm()$resp,topicsProb)
    k <- as.numeric(strsplit(as.character(input$topics1), "_")[[1]][2])
    text_temp<-resp%>%filter(topicsProb == k)
    resp_token<-sentenceParse(as.character(text_temp$text), docId = rownames(text_temp))
    Y<-resp_token %>% bind_lexrank(sentence,docId,sentenceId, threshold = input$s_threshold, level = 'sentences')
    A<-Y%>%
      group_by(docId)%>%
      summarize(Score=mean(lexrank),resp=paste(sentence, collapse=""))%>%
      arrange(desc(Score))%>%
      ungroup()%>%
      top_n(input$s_nbr,Score)

    print("dataframing the table")
    response <- data.frame(top_responses = A$resp)
  })
  
  # Suggest the best number of topics
  observeEvent(input$testButton, {
    output$elbow_curve <- renderPlotly({
      
      progress <- shiny::Progress$new()
      progress$set(message = "Computing Elbow Curve...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      
      if (is.null(input$comment)) {return()}
      if (is.null(input$gram_range[1])) {return()}
      if (is.null(input$gram_range[2])) {return()}
      text <- data()[, c(input$comment)]
      dtm <- prepare_dtm(text, TRUE, input$gram_range[1], input$gram_range[2], stopwords())$dtm
      f_dtm <- filter_dtm(dtm, 2)
      r_dtm <- remove_zero(f_dtm)$dtm
      optimal_k(r_dtm)
    })
  })
  
  # show top terms for a chosen topic
  output$top_terms <- renderTable({
    cat("printing top terms...")
    if (is.null(data())) {return()}
    as.data.frame(terms(result()$CTM, 10))
  })
  
  # tf_idf
  output$gg1 <- renderPlot({
    
    print("assigning topics")
    topics <- as.data.frame(topics(result()$CTM, 1))
    topics$id <- as.numeric(rownames(topics))
    colnames(topics) <- c('topic', 'id')
    
    print('reshaping dtm to tidy')
    tidy <- tidy(dtm()$dtm2) %>%
      mutate(id = as.numeric(document)) %>%
      inner_join(topics)
    
    print ("calculating tf_idf")
    tf_idf <- tidy %>%
      mutate(topic = as.factor(topic)) %>%
      select(topic, term) %>%
      group_by(topic) %>%
      count(term, sort = TRUE) %>%
      bind_tf_idf(term, topic, n) %>%
      arrange(desc(tf_idf)) %>%
      group_by(topic) %>%
      mutate(rank = row_number(tf_idf)) %>%
      filter(n >= 5) %>%
      ungroup()
    
    print("plotting things out")
    tf_idf %>%
      # filter(Csat %in% c(1, 10)) %>%
      group_by(topic) %>%
      top_n(12, rank) %>%
      mutate(word = reorder(term, tf_idf)) %>%
      ggplot(aes(word, tf_idf, fill = topic)) +
      geom_bar(alpha = 0.8, stat = "identity", show.legend = FALSE) +
      facet_wrap(~ topic, scales = "free", shrink = FALSE) +
      ylab("tf-idf") +
      coord_flip()
    
    
  })
  
  ################################ Metadata Information ##################################
  # show the metadata info
  output$filedf <- renderTable({
    if (is.null(data())) {return()}
    data.frame(summary(data()))
  })
  
  # data content  
  output$contents <- renderTable({
    if (is.null(data())) {return()}
    head(data())
  })
  
  # wordcloud
  observeEvent(input$goCloud, {
    output$wordcloud <- renderPlot({
      
      progress <- shiny::Progress$new()
      progress$set(message = "Plotting wordcloud...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      
      print("plotting wordcloud....")
      # m <- as.matrix(dtm()$dtm2)
      # data()$t = str_replace_all(data()[, input$comment],"[^[:graph:]]", " ") 
      # wordcloud(data()[data()[, input$group2] == input$group_options2,]$t)
      data <- data()
      data$t = str_replace_all(data[, input$comment],"[^[:graph:]]", " ")
      pal <- brewer.pal(6,"Dark2")
      pal <- pal[-(1)]
      data_sub <- data[data[, input$group2] == input$group_options2,]$t
      data_sub <- VCorpus(VectorSource(data_sub))
      data_sub <- tm_map(data_sub, removePunctuation)
      data_sub <- tm_map(data_sub, function(x)removeWords(x, c(stopwords(), tm::stopwords("english"), tm::stopwords('SMART'))))
      
      # wordcloud(data_sub,c(8,.3),2,100,TRUE,TRUE,.15,pal)
      wordcloud(data_sub, scale = c(8, .2), min.freq = 2, max.words = 100, colors = pal)
    })
  })

  
  ################################# Exploratory Analysis #################################
  tidy1 <- reactive({
    print(input$comment)
    
    dat <- data.frame(x = numeric(0), y = numeric(0))
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Shaping to tidy...", value = 0)

    # Number of times we'll go through the loop
    n <- 10

    for (i in 1:n) {
      # Each time through the loop, add another row of data. This is
      # a stand-in for a long-running computation.
      dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))

      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", i))

      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }
    
    if(is.null(data())) {return()}
    if(is.null(input$comment)) {return()}
    
    print("tidy the data into one gram tidy shape")
    tidy <- shape_to_tidy(data(), input$comment, input$gram)
  })
  
  
  
  tf_idf1 <- eventReactive(input$goTidy1, {
    
    dat <- data.frame(x = numeric(0), y = numeric(0))
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Calculating word importance...", value = 0)

    # Number of times we'll go through the loop
    n <- 10

    for (i in 1:n) {
      # Each time through the loop, add another row of data. This is
      # a stand-in for a long-running computation.
      dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))

      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", i))

      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }
    
    if (is.null(input$group1)) {return()}
    print(input$group1)
    print("calculating one gram tf_idf")
    tf_idf <- shape_to_tfidf(tidy1(), input$group1, input$threshold1)
    return(tf_idf)
  })
  
  # tf_idf2 <- eventReactive(input$goTidy2, {
  #   if (is.null(input$group2)) {return()}
  #   print(input$group2)
  #   print("calculating bigram tf_idf")
  #   tf_idf <- shape_to_tfidf(tidy2(), input$group1, input$threshold1)
  #   return(tf_idf)
  # })
  # 
  # group_count1 <- reactive({
  #   if (is.null(tf_idf1())) {return()}
  #   print("counting words for each group")
  #   tf_idf1() %>%
  #     group_by(group) %>%
  #     count(words = n()) %>%
  #     arrange(desc(words))
  # })
  
  # group_count2 <- reactive({
  #   if (is.null(tf_idf2())) {return()}
  #   print("counting words for each group")
  #   tf_idf2() %>%
  #     group_by(group) %>%
  #     count(words = n()) %>%
  #     arrange(desc(words))
  # })
  
  output$one_gram <- renderPlot({
    print("printing tf_idf")
    if (is.null(tf_idf1())) {return()}
    if (is.null(input$group_options)) {return()}
    # print()
    print(input$group_options)
    tf_idf1() %>%
      filter(group %in% input$group_options) %>%
      # filter(group %in% as.vector(as.matrix(head(group_count1(), input$howmany1)[,'group']))) %>%
      group_by(group) %>%
      top_n(12, rank) %>%
      mutate(word = reorder(word, tf_idf)) %>%
      ggplot(aes(word, tf_idf, fill = group)) +
      geom_bar(alpha = 0.8, stat = "identity", show.legend = FALSE) +
      facet_wrap(~ group, scales = "free", shrink = FALSE) +
      ylab("tf-idf") +
      coord_flip()
  })
  
  # output$bigram <- renderPlot({
  #   print("printing tf_idf")
  #   if (is.null(tf_idf2())) {return()}
  #   tf_idf2() %>%
  #     filter(group %in% as.numeric(input$group_options)) %>%
  #     # filter(group %in% as.vector(as.matrix(head(group_count2(), input$howmany2)[,'group']))) %>%
  #     group_by(group) %>%
  #     top_n(12, rank) %>%
  #     mutate(word = reorder(word, tf_idf)) %>%
  #     ggplot(aes(word, tf_idf, fill = group)) +
  #     geom_bar(alpha = 0.8, stat = "identity", show.legend = FALSE) +
  #     facet_wrap(~ group, scales = "free", shrink = FALSE) +
  #     ylab("tf-idf") +
  #     coord_flip()
  # })

  ################################# UI Setup ##############################################
  # UI that enable user to select a topic to plot
  output$topic_choices <- renderUI({
    n = input$nbr_topics
    selectInput('topics', 'choose a topic to plot', 
                paste0('topic_', seq_along(1:n)))
  })
  
  output$topic_choices1 <- renderUI({
    n = input$nbr_topics
    selectInput('topics1', 'choose a topic to plot', 
                paste0('topic_', seq_along(1:n)))
  })
  
  output$gram <- renderUI({
    selectInput('gram', 
                'choose one or bigram', 
                c(1, 2))
  })
  
  # UI that indicates which column is the Comment/Openend/Text
  output$Comment <- renderUI({
    columns <- colnames(data())
    selectInput('comment',
                'Choose text column', 
                choices = columns)
  })
  
  output$s_threshold <- renderUI({
    sliderInput("s_threshold", "Page Rank Threshold",
                min = 0.1, max = 1, value = 0.15)
  })
  
  output$s_nbr <- renderUI({
    sliderInput("s_nbr", "Number of sentence to show",
                min = 1, max = 10, value = 5)
  })
  
  output$group1 <- renderUI({
    columns <- colnames(data())
    selectInput('group1', 
                'Group data into',
                choices = columns)
  })
  
  output$group2 <- renderUI({
    columns <- colnames(data())
    selectInput('group2', 
                'Group data into',
                choices = columns)
  })
  
  output$group3 <- renderUI({
    columns <- colnames(data())
    selectInput('group3', 
                'Group data into',
                choices = columns)
  })
  
  output$threshold3 <- renderUI({
    sliderInput("threshold3", "Word Frequency Filter",
                min = 10, max = 200, value = 50)
  })
  
  output$threshold1 <- renderUI({
    sliderInput("threshold1", "Word Frequency Filter",
                min = 1, max = 30, value = 1)
  })
  
  output$threshold2 <- renderUI({
    sliderInput("threshold2", "Word Frequency Filter",
                min = 1, max = 30, value = 1)
  })
  
  output$show_how_many1 <- renderUI({
    sliderInput("howmany1", "Please select how many plots you want to see:",
                min = 1, max = 6, value = 2)
  })
  
  output$show_how_many2 <- renderUI({
    sliderInput("howmany2", "Please select how many plots you want to see:",
                min = 1, max = 6, value = 2)
  })
  
  output$min_gram <- renderUI({
    sliderInput('gram_range', 
                'Gram range', 
                min = 1, max = 6, value = c(1, 2))
    # selectInput()
  })
  
  
  output$net_threshold <- renderUI({
    sliderInput('net_threshold',
                'Minimum frequency of words',
                min = 1, max = 50, value = 8)
  })
  
  
  output$group_options <- renderUI({
    # unique_values <- as.character(unique(data()[, input$group1])[order(unique(data()[,input$group1]))])
    if (is.null(input$group1)) {return()}
    print(paste("group1 = ", input$group1))
    unique_values <- unique(data()[, input$group1])
    checkboxGroupInput("group_options", label = h3("group options"), 
                       choices = unique_values, selected = unique_values[1])
  })
  
  output$group_options2 <- renderUI({
    # unique_values <- as.character(unique(data()[, input$group2])[order(unique(data()[,input$group2]))])
    unique_values <- unique(data()[, input$group2])
    selectInput("group_options2", label = h3("group options"), 
                       choices = unique_values, selected = unique_values[1])
  })
  
  output$group_options3 <- renderUI({
    # unique_values <- as.character(unique(data()[, input$group3])[order(unique(data()[,input$group3]))])
    if (is.null(input$group3)) {return()}
    unique_values <- unique(data()[, input$group3])
    selectInput("group_options3", label = h3("group options"), 
                choices = unique_values, selected = unique_values[1])
  })

  observeEvent(input$goButton, {
    if (is.null(topic())) {return()}
    output$network_proxy_nodes <- renderVisNetwork({
      cat("drawing the plot...")
      dat <- data.frame(x = numeric(0), y = numeric(0))
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Making plot", value = 0)
      
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      
      dataset <- topic()
      nodes <- create_node(dataset$topic, dataset$prob, input$net_threshold)
      edges <- create_edge(dataset$topic, nodes)
      nodes_clean <- nodes %>%
        filter(nodes$label %in% edges$from)
      lnodes <- data.frame(label = c("low probability", "high probability"),
                           shape = c( "ellipse"), color = c("red", "pink"),
                           id = 1:2)
  
      visNetwork(nodes_clean, edges, main = "Words Occurence Network", height = "1200px", width = "100%") %>%
        visOptions(highlightNearest = TRUE,
                   nodesIdSelection = TRUE) %>%
        visNodes(font = '35px arial black') %>%
        visLegend(addNodes = lnodes, useGroup = FALSE, main = 'Legend') %>%
        visInteraction(navigationButtons = TRUE)%>%
        visPhysics(stabilization = TRUE) %>%
        # visIgraphLayout() %>%
        visOptions(manipulation = TRUE, highlightNearest = TRUE)
  
    })
  })
  
  tidy_edge <- reactive({
    print("making edge")
    if(is.null(input$group3)) {return()}
    if(is.null(input$group_options3)) {return()}
    if(is.null(input$comment)) {return()}
    tidy2 <- shape_to_tidy(data(), input$comment, 2)
    edge <- tidy2 %>%
      separate(word, c("word1", "word2"), sep = " ") %>%
      mutate_(group = interp(quote(as.factor(var)), var = as.name(input$group3))) %>%
      group_by(group) %>%
      count(word1, word2, sort = TRUE) %>%
      mutate(from = word1) %>%
      mutate(to = word2) %>%
      mutate(width = log(n)) %>%
      filter(group %in% input$group_options3)
      # filter(n > 20) %>%
    # print("select group for edge")
    # edge[edge[, input$group3] == input$group_options3,]
  })

  
  tidy_node <- reactive({
    print("making node")
    if(is.null(input$group3)) {return()}
    if(is.null(input$group_options3)) {return()}
    if(is.null(input$comment)) {return()}
    tidy1 <- shape_to_tidy(data(), input$comment, 1)
    node <- tidy1 %>%
      mutate_(group = interp(quote(as.factor(var)), var = as.name(input$group3))) %>%
      group_by(group) %>%
      count(word, sort = TRUE) %>%
      mutate(label = word) %>%
      mutate(id = word) %>%
      mutate(value = log(n)) %>%
      filter(n > input$threshold3) %>%
      filter(group %in% input$group_options3)
    # print("select group for node")
    # node1 <- node[node[, input$group3] == input$group_options3,]
    print("clean the node")
    clean_node <- node %>%
      filter(node$label %in% tidy_edge()$from)
  })
  
  observeEvent(input$goAgain, {
    output$network1 <- renderVisNetwork({
      cat("drawing the plot...")
      dat <- data.frame(x = numeric(0), y = numeric(0))
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Making plot", value = 0)
      
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      
      visNetwork(tidy_node(), tidy_edge(), main = "Words Occurence Network", height = "1200px", width = "100%") %>%
        visOptions(highlightNearest = TRUE,
                   nodesIdSelection = TRUE) %>%
        visNodes(font = '35px arial black') %>%
        # visLegend(addNodes = lnodes, useGroup = FALSE, main = 'Legend') %>%
        visInteraction(navigationButtons = TRUE)%>%
        visPhysics(stabilization = TRUE) %>%
        # visIgraphLayout() %>%
        visOptions(manipulation = TRUE, highlightNearest = TRUE)
      
        
      })
    })
  
  observeEvent(input$goButton, {
    if (is.null(topic())) {return()}
    output$network <- renderVisNetwork({
      cat("drawing the plot...")
      dat <- data.frame(x = numeric(0), y = numeric(0))
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Making plot", value = 0)
      
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      
      dataset <- topic()
      nodes <- create_node(dataset$topic, dataset$prob, input$net_threshold)
      edges <- create_edge(dataset$topic, nodes)
      nodes_clean <- nodes %>%
        filter(nodes$label %in% edges$from)
      lnodes <- data.frame(label = c("low probability", "high probability"),
                           shape = c( "ellipse"), color = c("red", "pink"),
                           id = 1:2)
      
      visNetwork(nodes_clean, edges, main = "Words Occurence Network", height = "1200px", width = "100%") %>%
        visOptions(highlightNearest = TRUE,
                   nodesIdSelection = TRUE) %>%
        visNodes(font = '35px arial black') %>%
        visLegend(addNodes = lnodes, useGroup = FALSE, main = 'Legend') %>%
        visInteraction(navigationButtons = TRUE)%>%
        visPhysics(stabilization = TRUE) %>%
        # visIgraphLayout() %>%
        visOptions(manipulation = TRUE, highlightNearest = TRUE)
      
    })
  })
}



ui <- tagList(
  navbarPage(
    theme = "cerulean",  # <--- To use a theme, uncomment this
    "Text Mining Tool",
    tabPanel("Data Upload and Summary",
             sidebarPanel(
               fileInput('file1', 'Upload your file',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv')),
               radioButtons('sep1', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ','),
               uiOutput('Comment'),
               fileInput('file2', 'Upload your stopwords',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv'))
               # radioButtons('quote', 'Quote',
               #              c(None='',
               #                'Double Quote'='"',
               #                'Single Quote'="'"),
               #              '"')
             ),
             mainPanel(
               if(is.null(data())) 
                 h5("please upload your dataset") 
               else
                 tabsetPanel(tabPanel("About file", tableOutput("filedf")), 
                             tabPanel("Sample data", tableOutput('contents')))
                             # tabPanel("topics", tableOutput('top_terms')))
             )
    ),
    tabPanel("Exploratory Analysis", 
             tabsetPanel(
               tabPanel('Important Words',
                        sidebarPanel(
                          uiOutput('group1'),
                          uiOutput('gram'),
                          uiOutput('group_options'),
                          uiOutput('threshold1'),
                          actionButton('goTidy1', 'Run TF_IDF')
                        ),
                        mainPanel(plotOutput('one_gram'))),
               tabPanel('wordcloud', 
                        sidebarPanel(
                          uiOutput('group2'),
                          uiOutput('group_options2'),
                          # uiOutput('show_how_many2'),
                          actionButton('goCloud', 'Draw the wordcloud')
                        ),
                        mainPanel(plotOutput('wordcloud'))),
               tabPanel('Network View',
                        sidebarPanel(
                          uiOutput('group3'),
                          uiOutput('group_options3'),
                          uiOutput('threshold3'),
                          actionButton('goAgain', 'Draw the network')
                        ),
                        mainPanel(
                          visNetworkOutput("network1")
                        )
                      
               )
             )
    ),
    tabPanel("Topic Modeling",
             sidebarPanel(
               selectInput('nbr_topics',
                           'please select number of topics you want to have',
                           2:20),
               uiOutput('min_gram'),
               actionButton("goButton", "Run Topic Modeling!"),
               downloadButton('downloadData', 'Download')
               # uiOutput('topic_choices'),
               # uiOutput('net_threshold')
               # uiOutput('s_threshold'),
               # uiOutput('s_nbr'),
               # actionButton("testButton", "Click to see the optimal number of topics")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel('Term Network Analysis', 
                          fluidRow(
                            column(width = 6, uiOutput('topic_choices')),
                            column(width = 6, uiOutput('net_threshold'))
                            ),
                          fluidRow(width = 12, 
                                   visNetworkOutput("network_proxy_nodes", height = "400px"))
                          ),
                 tabPanel('Important Sentences', 
                          fluidRow(
                            column(width = 4, uiOutput('topic_choices1')),
                            column(width = 4, uiOutput('s_threshold')),
                            column(width = 4, uiOutput('s_nbr'))
                            ),
                          fluidRow(width = 12, tableOutput("top_sentence"))
                          ),
                 tabPanel('Important Terms', plotOutput('gg1')),
                 tabPanel('Topic Distribution', plotlyOutput('pie_chart')),
                 tabPanel('Elbow Curve', 
                          fluidRow(
                            column(width = 6, 
                                   actionButton("testButton", "Click to see the optimal number of topics"))
                            ),
                          fluidRow(
                            width = 12, plotlyOutput('elbow_curve'))
                 )
               )

             )
    )
  )
)

shinyApp(ui = ui, server = server)

