library(M4metalearning)

calc_errors_new <- function(dataset) {
  
  total_snaive_errors <- c(0,0)
  for (i in 1:length(dataset)) {
    tryCatch({
      lentry <- dataset[[i]]
      insample <- lentry$x
      
      #extrac forecasts and attach the snaive for completion
      ff <- lentry$ff
      ff <- rbind(ff, snaive_forec(insample, lentry$h))
      
      frq <- frq <- stats::frequency(insample)
      insample <- as.numeric(insample)
      outsample <- as.numeric(lentry$xx)
      masep <- mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))
      
      
      repoutsample <- matrix(
        rep(outsample, each=nrow(ff)),
        nrow=nrow(ff))
      ## |A_t-F_t|/(|A_t|+|F_t|)/2*100
      smape_err <- 200*abs(ff - repoutsample) / (abs(ff) + abs(repoutsample))
      ## |A_t-F_t|/masep
      mase_err <- abs(ff - repoutsample) / masep
      ## |A_t-F_t|/|A_t|
      mape_err <- abs(ff- repoutsample)/abs(repoutsample)
      ## |A_t-F_t|
      mae_err <- abs(ff- repoutsample) ##/as.numeric(c(1,1,1,1))
      
      
      lentry$snaive_mase <- mase_err[nrow(mase_err), ]
      lentry$snaive_smape <- smape_err[nrow(smape_err),]
      
      lentry$mase_err <- mase_err[-nrow(mase_err),]
      lentry$smape_err <- smape_err[-nrow(smape_err),]
      
      lentry$mape_err <- mape_err[-nrow(mape_err),]
      lentry$mae_err <- mae_err[-nrow(mae_err),]
      dataset[[i]] <- lentry
      total_snaive_errors <- total_snaive_errors + c(mean(lentry$snaive_mase),
                                                     mean(lentry$snaive_smape))
    } , error = function (e) {
      print(paste("Error when processing OWIs in series: ", i))
      print(e)
      e
    })
  }
  total_snaive_errors = total_snaive_errors / length(dataset)
  avg_snaive_errors <- list(avg_mase=total_snaive_errors[1],
                            avg_smape=total_snaive_errors[2])
  
  
  for (i in 1:length(dataset)) {
    lentry <- dataset[[i]]
    #dataset[[i]]$errors <- 0.5*(rowMeans(lentry$mase_err)/avg_snaive_errors$avg_mase +
    #                              rowMeans(lentry$smape_err)/avg_snaive_errors$avg_smape)
    dataset[[i]]$errors <- rowMeans(lentry$mae_err^2)
  }
  attr(dataset, "avg_snaive_errors") <- avg_snaive_errors
  dataset
}

THA_features_new <-function(dataset, n.cores=1) {
  list_process_fun <- lapply
  cl = -1
  require(tsfeatures)
  dataset_feat <- list_process_fun(dataset,
                                 function (serdat) {
                                   tryCatch({
                                     #additional features from Talagala, Hyndman, Athanasopoulos 2018
                                     featrow <-
                                       tsfeatures::tsfeatures(
                                         serdat$x,
                                         features = c(
                                           "acf_features",
                                           "arch_stat",
                                           "crossing_points",
                                           "entropy",
                                           "flat_spots",
                                           "heterogeneity_tsfeat_workaround",
                                           "holt_parameters",
                                           "hurst",
                                           "lumpiness",
                                           "nonlinearity",
                                           "pacf_features",
                                           "stl_features",
                                           "stability",
                                           #"hw_parameters_tsfeat_workaround",
                                           "unitroot_kpss",
                                           "unitroot_pp"
                                         )
                                       )
                                     
                                     
                                     #additional features
                                     series_length <- length(serdat$x)
                                     
                                     featrow <- tibble::add_column(
                                       featrow,
                                       "series_length" = series_length)
                                     
                                     featrow[is.na(featrow)] <-
                                       0 #SET NAs TO 0 ?
                                     
                                     
                                     serdat$features <- featrow
                                     serdat
                                   }, error = function(e) {
                                     print(e)
                                     return(e)
                                   })
                                 })

dataset_feat
}

true_data <- read.table("~/Covid-Forecast/IHME/Rhode Island_true_data.txt", quote="\"", comment.char="")
truex = as.list(true_data)
x <- ts(truex$V1, frequency = 1)
ele <- list(list(x=x[1:67], h = 4)) ##0524

for(i in 0:13){
  #print(67+i)
  true_data <- read.table("~/Covid-Forecast/IHME/Rhode Island_true_data.txt", quote="\"", comment.char="")
  truex = as.list(true_data)
  temp <- ts(truex$V1, frequency = 1)
  #print(temp[1:(67+i)])
  ele <- append(ele, list(list(x=temp[1:(69+i)], h = 4))) ##0607
}
#ele <- append(ele, list(list(x=x[1:67], h = 4))) ##0607
#ele <- append(ele, list(list(x=x[1:68], h = 4))) ##0614
#ele <- append(ele, list(list(x=x[1:69], h = 4))) ##0621
#ele <- append(ele, list(list(x=x[1:70], h = 4))) ##0628
#ele <- append(ele, list(list(x=x[1:71], h = 4))) ##0705
#ele <- append(ele, list(list(x=x[1:72], h = 4))) ##0712
#ele <- append(ele, list(list(x=x[1:73], h = 4))) ##0719
#ele <- append(ele, list(list(x=x[1:74], h = 4))) ##0726
#ele <- append(ele, list(list(x=x[1:75], h = 4))) ##0802
#ele <- append(ele, list(list(x=x[1:76], h = 4))) ##0809
#ele <- append(ele, list(list(x=x[1:77], h = 4))) ##0816

#state_names <- as.list(c())
#state_names <- as.list(c('Kentucky','Missouri','Texas','New Mexico'))
#state_names <- as.list(c("Illinois","Georgia","Oregon","Florida","New Hampshire","Minnesota",'Ohio','Nevada',"Arizona",'Texas','Kentucky','Missouri','Rhode Island','New Mexico'))
#state_names <- as.list(c('South Dakota','North Dakota'))
state_names <- as.list(c('Connecticut', 'Maryland', 'Nevada', 'New Hampshire', 'New Jersey', 'Pennsylvania'))
#state_names <- as.list(c('Kansas', 'Tennessee', 'West Virginia', 'Montana'))
#state_names <- as.list(c('Pennsylvania','Arkansas','Wisconsin','Idaho','North Carolina', 'Indiana','Iowa'))
#state_names <- as.list(c("Virginia", 'Mississippi','Alabama','South Carolina','Nebraska','Oklahoma'))
#state_names <- as.list(c("Connecticut",'Michigan','New Jersey','Colorado'))
#state_names <- as.list(c("Louisiana", 'Delaware','Massachusetts'))
#state_names <- as.list(c("Illinois","Georgia","Oregon","Rhode Island","New Hampshire","Alabama","Delaware","Maryland","United States of America"))
#state_names <- as.list(c("South Carolina","Ohio","Nevada","Missouri","New Mexico","Oklahoma","Utah","Kentucky","Arizona", "Minnesota"))
#state_names <- as.list(c("Maryland","Delaware","Alabama","Illinois","New Hampshire","Georgia","Iowa","Oregon","Nebraska"))
#state_names <- as.list(c('Texas','Washington','United States of America',"Florida","New York",'Alabama', 'Arizona', 'Arkansas', 'Colorado', 'Connecticut','Delaware','Georgia','Idaho', 'Illinois','Indiana', 'Iowa',
#                         'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', "Massachusetts", 'Michigan','Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico',  'North Carolina', 'North Dakota',
#                         'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee','Utah', 'Vermont', 'Virginia',  'West Virginia', 'Wisconsin', 'Wyoming'))
#for(state in levels(state_names$V1)){
for(state in state_names){
  #if (state==testing_state){
  ## import the full lists
  state_filename <- paste("~/Covid-Forecast/IHME/", paste(state,"_true_data.txt", sep = ""), sep="")
  #print(state_filename)
  true_data <- read.table(state_filename, quote="\"", comment.char="")
  truex = as.list(true_data)
  ## extract time series from full lists
  
  for(i in 0:15){
    
    #if(i==0){
    #temp1 <- ts(truex$V1, frequency = 1)
    #ele<-list(list(x=temp1[1:65], h=4))
    #}
    if(i!=1){
      temp2 <- ts(truex$V1, frequency = 1)
      xstart <- sample(1:10, 1)
      ele <- append(ele, list(list(x=temp2[1:(67+i)], h=4))) ##target date+4
    }
  }
  #}
}
ele <- temp_holdout(ele)

state_names <- c('Rhode Island', state_names)

## fill in state forecasts
i = 1
date_list = list(c("05-24","06-07","06-14","06-21","06-28","07-05","07-12","07-19","07-26","08-02","08-09","08-16","08-23","08-30","09-06"))
for(state in state_names){
  #if (state==testing_state){
  #j = 1
  for(d in date_list){
    state_forecast_filename <- paste(paste(paste("~/Covid-Forecast/", paste(state,"_forecast_list_",sep = ""), sep=""),d, sep = ""),".txt",sep = "")
    
    for(k in 1:15){
      #print(state_forecast_filename[k])
      temp <- as.matrix(read.table(state_forecast_filename[k], sep = ","))
      #print(temp)
      ele[[i]]$ff <- temp
      i = i+1
      #print(i)
    }
  }
  #}
}

ele <- calc_errors(ele)
ele<- THA_features(ele)
train_data_covid <- create_feat_classif_problem(ele)
set.seed(1345)
meta_model_covid <- train_selection_ensemble(train_data_covid$data, train_data_covid$errors)

## visualize and validate xgboost model
mat <- xgboost::xgb.importance (feature_names = colnames(train_data_covid$data),model = meta_model_covid)
xgboost::xgb.plot.importance(importance_matrix = mat[1:20], cex=1.0)

############################################### testing #################################################################
true_data <- read.table("~/Covid-Forecast/IHME/Pennsylvania_true_data.txt", quote="\"", comment.char="")
truex = as.list(true_data)

x <- ts(truex$V1, frequency = 1)
test <- list(list(x=x[1:83], h=4)) ##0913
test <- append(test, list(list(x=x[1:84], h = 4))) ##0920
test <- append(test, list(list(x=x[1:85], h = 4))) ##0927
test <- append(test, list(list(x=x[1:86], h = 4))) ##1004
test <- append(test, list(list(x=x[1:87], h = 4))) ##1011
test <- temp_holdout(test)



date_list = c("09-13","09-20","09-27","10-04","10-11")
i=1
for(d in date_list){
  state_forecast_filename <- paste(paste("~/Covid-Forecast/Pennsylvania_forecast_list_",d, sep = ""),".txt",sep = "")
  #print(state_forecast_filename)
  temp <- as.matrix(read.table(state_forecast_filename, sep = ","))
  test[[i]]$ff <- temp
  i = i+1
}

test <- THA_features(test, n.cores=1)
test_data <- create_feat_classif_problem(test)
preds <- predict_selection_ensemble(meta_model_covid, test_data$data)
head(preds)
test <- ensemble_forecast(preds, test)


#IHME<- c(3.469660735050057, 3.776959051242368, 4.0331975183276905, 4.2496625594276045, 4.415196302306807, 4.5770007729984306, 4.759997541510252, 4.974343219549579, 5.160884580501815, 5.3594575008650285, 5.494961422159894, 5.5141538692375605, 5.603284303244241, 5.6987374491095375) ##Arizona
#IHME <- c(13316, 14863, 15592, 15704, 15059, 14503, 13930, 13311, 12692, 11781, 10897)/331893745 ##0823 national

#IHME <- c(2.0978316490693776, 2.8680456593956705, 3.6173964540777694, 4.312210771536068, 4.83547478049262, 5.172225616422, 5.3077384395192375, 5.311749194825172, 5.164862615025899, 4.8251533980160985, 4.39705689281858, 3.9480807609347934, 3.5350322767723537, 3.2308162671417855, 3.060216506872221) ## WA
IHME <- c(1.8272752860274966, 2.408560706558843, 2.934336643845496, 3.33815809821928, 3.574438712675184, 3.6585498751766026, 3.6191919452438706, 3.5149173364198774, 3.411132306073148, 3.180875469279192, 2.9736966126266493, 2.811292141595169, 2.6876069312283177, 2.599751168744372, 2.5469702253712465) ##CA



date_list = c("09-13","09-20","09-27","10-04","10-11")



for(i in 1:5){
  y_avg <- colMeans(test[[i]]$ff)
  
  test[[i]]$ff <- rbind(test[[i]]$ff, y_avg)
  
  test[[i]]$ff <- rbind(test[[i]]$ff, IHME[i:(i+3)])
  filename = paste(paste(paste("~/Covid-Forecast/", "Pennsylvania_CDC_Ensemble_list_", sep=""),date_list[i], sep = ""),".txt",sep = "")
  tempcdc=as.matrix(read.table(filename,sep=","))
  CDC_ensemble = tempcdc[1,]
  test[[i]]$ff <- rbind(test[[i]]$ff, CDC_ensemble)
  
  test[[i]]$ff <- rbind(test[[i]]$ff, test[[i]]$y_hat)
}
result <-calc_errors_new(test)

mse_avg = rep(0, 9)
for (i in 1:5){
  filename = paste(i, "US_smape.csv", sep = "")
  write.csv(result[[i]]$mase_err, filename)
  
  mse_avg = mse_avg+result[[i]]$errors
}
mse_avg = sqrt(mse_avg/5)
print(mse_avg)
