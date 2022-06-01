library(M4metalearning)

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



true_data <- read.table("~/Covid-Forecast/IHME/United States of America_true_data.txt", quote="\"", comment.char="")
truex = as.list(true_data)
## need to change the index here
x <- ts(truex$V1, frequency = 1)
ele <- list(list(x=x[1:61], h = 4)) ##0524

state_names <- as.list(c("California","Texas","New York","Washington",'Alabama', 'Arizona', 'Arkansas', 'Colorado', 'Connecticut','Delaware', 'Florida','Georgia','Idaho', 'Illinois','Indiana', 'Iowa',
                         'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan','Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico',  'North Carolina', 'North Dakota',
                         'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee','Utah', 'Vermont', 'Virginia',  'West Virginia', 'Wisconsin', 'Wyoming','District of Columbia', 'Alaska'))

for(state in state_names){
  #if (state==testing_state){
  ## import the full lists
  state_filename <- paste("~/Covid-Forecast/IHME/", paste(state,"_true_data.txt", sep = ""), sep="")
  #print(state_filename)
  true_data <- read.table(state_filename, quote="\"", comment.char="")
  truex = as.list(true_data)
  temp2 <- ts(truex$V1, frequency = 1)
  ## extract time series from full lists
  ele <- append(ele, list(list(x=temp2[1:61], h=4))) ##target date+4
}
ele <- temp_holdout(ele)
ele<- THA_features_new(ele)

ele<-calc_forecasts(ele, forec_methods())
ele <- calc_errors(ele)

train_data_covid <- create_feat_classif_problem(ele)

X = train_data_covid$data
write.table(X, file = 'covid_features.txt')

X.pca <- prcomp(X, center = TRUE,rank=2)
pcs <- X.pca$x
rownames(pcs)<-state_names
model=kmeans(X,3)
library(cluster)
clusplot(pcs,model$cluster, labels = 2)
