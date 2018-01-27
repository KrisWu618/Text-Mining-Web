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
  df <- data.frame(nbr_topics = sequ, log_likelihood = hm_many)
  plot_ly(df, x = ~nbr_topics, y = ~log_likelihood) %>%
    add_lines()
   
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
