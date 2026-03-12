setwd("C:/Users/simon/OneDrive/Desktop/licenta_cod")

library(clusterCrit)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
library(ggplot2)
library(plotly)
library(dplyr)
library(countrycode)
library(ggthemes)
library(moments)
library(corrplot)
library(factoextra)
library(FactoMineR)
library(psych)
library(nFactors)
library(ca)
library(cluster)
library(dendextend) 
library(Hmisc)
library(car)
library(NbClust)
library(plotly)
library(pls)
library(boot)     
library(stats) 
library(lmtest)
library(ResourceSelection)
library(pROC)
library(MASS)

##citirea datelor
data<-read.table(file="date5.txt", header=TRUE, dec=".", row.names=1)
dim(data)

##eliminare outliers
detect_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  x < (q1 - 3.5 * iqr) | x > (q3 + 3.5 * iqr)
}
outliers_matrix <- sapply(data, detect_outliers)
rows_out <- which(apply(outliers_matrix, 1, any))
print(rownames(data)[rows_out]) # "DK" "DE" "EE" "LU"
date <- data[-rows_out, ]
data1 <-data[-rows_out, ]
corelatie<-cor(date)
corrplot::corrplot(corelatie,method='number', type="upper")
date_bun<-date[,-c(2,5,8)]

corelatie1<-cor(date_bun)
corrplot::corrplot(corelatie1,method='number', type="upper")

#Statistici descriptive
summary(date_bun)
skewness(date_bun)
kurtosis(date_bun)

#boxplot-uri
par(mfrow=c(3,3))
boxplot(date_bun$X1, main = "Boxplot pentru variabila X1", col = "skyblue")
boxplot(date_bun$X3, main = "Boxplot pentru variabila X3", col = "skyblue")
boxplot(date_bun$X4, main = "Boxplot pentru variabila X4", col = "skyblue")
boxplot(date_bun$X6, main = "Boxplot pentru variabila X6", col = "skyblue")
boxplot(date_bun$X7, main = "Boxplot pentru variabila X7", col = "skyblue")
boxplot(date_bun$X9, main = "Boxplot pentru variabila X9", col = "skyblue")
boxplot(date_bun$X10, main = "Boxplot pentru variabila X10", col = "skyblue")
boxplot(date_bun$X11, main = "Boxplot pentru variabila X11", col = "skyblue")
par(mfrow=c(1,1))

#matricea de corelație
rez <- rcorr(as.matrix(date_bun))
rez
corelatii_filtrate <- rez$r
corelatii_filtrate[rez$P > 0.05] <- NA 

#corelații semnificative
corrplot::corrplot(corelatii_filtrate, 
                   method = "color", 
                   type = "upper", 
                   na.label = " ", 
                   addCoef.col = "black", 
                   title = "Corelații semnificative (p < 0.05)")
#-------------------------------------------------------------------------------------------------
#ANALIZA COMPONENTELOR PRINCIPALE
#------------------------------------------------------------------------------------------------
#Standardizarea datelor
date_std<- scale(date_bun)
cor(date_std)
cov(date_std)

#Aplicarea ACP cu funcția princomp
acp <- princomp(date_std, cor = TRUE)
summary(acp)
fviz_eig(acp)
round(acp$loadings[, 1:3],2)

biplot(acp, main = "Biplot ACP")
fviz_pca_ind(acp)
fviz_pca_ind(acp,
             label     = "all",   
             repel     = TRUE,    
             labelsize = 6)     
fviz_pca_var(acp, col.var="contrib")
#scorurile acp
scoruri <-acp$scores
scoruri
cor(date_std, scoruri)

##ACP cu funcția PCA
acp2<-PCA(date_std,scale.unit=T, graph=F)
acp2
summary(acp2)
fviz_pca_var(acp2, col.var="contrib")
fviz_pca_ind(acp2) 
par(mfrow=c(1,3))
fviz_contrib(acp2, choice = "var", axes = 1, top = 10)  # pentru PC1
fviz_contrib(acp2, choice = "var", axes = 2, top = 10)  # pentru PC2
fviz_contrib(acp2, choice = "var", axes = 3, top = 10)  # pentru PC3

#scorurile acp2
pca_scores <- acp2$ind$coord
pca_scores

#grafic 3D al scorurilor acp
plot_ly(
  x = ~pca_scores[,1],
  y = ~pca_scores[,2],
  z = ~pca_scores[,3],
  type = "scatter3d",
  mode = "markers+text",
  text = rownames(date_bun), 
  marker = list(
    size = 5,
    color = ~pca_scores[,1], 
    colorscale = "Viridis",  
    colorbar = list(title = "Dim 1")
  ),
  textposition = "top center"
) %>%
  layout(scene = list(
    xaxis = list(title = "Dim 1"),
    yaxis = list(title = "Dim 2"),
    zaxis = list(title = "Dim 3")
  ))

#HARTA
scoruri_harta <- data.frame(acp$scores[, 1:3]) 
colnames(scoruri_harta) <- c("CP1", "CP2", "CP3")

scoruri_harta$iso2c <- rownames(date_bun)  
scoruri_harta$iso_a3 <- countrycode(scoruri_harta$iso2c, origin = "iso2c", destination = "iso3c")
scoruri_harta$iso2c <- trimws(rownames(date_bun))
scoruri_harta$iso_a3 <- countrycode(scoruri_harta$iso2c, origin = "iso2c", destination = "iso3c")
scoruri_harta %>% filter(iso2c == "FR") 
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- world %>% filter(continent == "Europe")
europe <- world %>% filter(iso_a3 %in% scoruri_harta$iso_a3)
map_data <- left_join(europe, scoruri_harta, by = "iso_a3")
map_data <- europe %>% 
  right_join(scoruri_harta, by = "iso_a3") %>%
  filter(!is.na(CP1))
map_data <- map_data %>% filter( !is.na(CP1))

h1 <- ggplot(map_data) +
  geom_sf(aes(fill = CP1)) +
  scale_fill_gradient(low = "#ffe5e5", high = "darkred", name = "Scor CP1") +
  ggtitle("CP 1 (Performanță economică și digitală)") +
  coord_sf(xlim = c(-10, 35), ylim = c(35, 72)) +
  theme_minimal(base_size = 10)

h2 <- ggplot(map_data) +
  geom_sf(aes(fill = CP2)) +
  scale_fill_gradient(low = "#eaffea", high = "darkgreen", name = "Scor CP2") +
  ggtitle("CP 2 (Presiune fiscală cu eficiență scăzută a inovației )") +
  coord_sf(xlim = c(-10, 35), ylim = c(35, 72)) +
  theme_minimal(base_size = 9)

h3 <- ggplot(map_data) +
  geom_sf(aes(fill = CP3)) +
  scale_fill_gradient(low = "#e0f0ff", high = "darkblue", name = "Scor CP3") +
  ggtitle("CP 3 (Antreprenoriat cantitativ, dar nu și inovativ) ") +
  coord_sf(xlim = c(-10, 35), ylim = c(35, 72)) +
  theme_minimal(base_size = 9)  
h1
h2
h3
#pregătesc scorurile acp pentru analiza cluster
scoruri <- as.data.frame(acp$scores)
cluster_data <- scoruri[, 1:3]

#-------------------------------------------------------------------------------------------------
#ANALIZA CLUSTER PE COMPONENTELE PRINCIPALE
#------------------------------------------------------------------------------------------------
d <- dist(scoruri)
d
ierarhie_ward <- hclust(d, method = "ward.D2")
plot(ierarhie_ward, main = "Dendogramă - Ward")
ierarhie_ward$height 

nr <- NbClust(scoruri, distance="euclidean", min.nc=2, max.nc=7, method="ward.D2", index="all")
nr$All.index
nr$Best.nc
table(nr$Best.n[1, ])

#varianta cu 2 clustere
rect.hclust(ierarhie_ward, k = 2, border = 2:5)
clase_ward <- cutree(ierarhie_ward, k = 2)
s_ward <- silhouette(clase_ward, d)
plot(s_ward, main = "Silhouette - Ward") #0,34

#varianta cu 3 clustere
rect.hclust(ierarhie_ward, k = 3, border = 2:5)
clase_ward1 <- cutree(ierarhie_ward, k = 3)
s_ward1 <- silhouette(clase_ward1, d)
plot(s_ward1, main = "Silhouette - Ward") #0,21
par(mfrow=c(1,1))

fviz_cluster(list(data = d, cluster = clase_ward))+
  labs(title = "Reprezentare vizuală a clusterelor prin metoda Ward")

#Elbow method
fviz_nbclust(cluster_data, kmeans, method = "wss") +
  labs(title = "Metoda Elbow")

#Silhouette method
fviz_nbclust(cluster_data, kmeans, method = "silhouette") +
  labs(title = "Metoda Silhouette")

# Analiza cluster cu algoritmul K-means
set.seed(123) 
kmeans_result1 <- kmeans(cluster_data, centers = 2, nstart = 25) 
cluster_data$cluster <- factor(kmeans_result1$cluster)
clase_km <- kmeans_result1$cluster
table(clase_km)
d_km <- dist(cluster_data)
fviz_cluster(kmeans_result1, d_km) +
  labs(title = "Reprezentarea vizuală a clusterelor prin metoda K-means")
s <- silhouette(clase_km,d_km)
plot(s)

set.seed(12)  
kmeans_result <- kmeans(cluster_data, centers = 3, nstart = 25) 
cluster_data$cluster <- factor(kmeans_result$cluster)
clase <- kmeans_result$cluster 
table(clase)
d_km1 <- dist(cluster_data)
fviz_cluster(kmeans_result, d_km1) +
  labs(title = "Reprezentarea vizuală a clusterelor prin metoda K-means")
s <- silhouette(clase,d_km1)
plot(s)
kmeans_result1$cluster
kmeans_result1$centers
kmeans_result1$totss 
kmeans_result1$withinss 
kmeans_result1$tot.withinss 
kmeans_result1$betweenss 

#Harta
cluster_data_harta <- data.frame(
  iso_a2 = c("BE", "FI", "AT", "NL","DK","SE","RO", "BG", "LV", "GR", "HU", "PL", "LT", "HR", "PT", "CZ", "SI", "IT", "ES", "CY", "FR"),
  cluster = c(rep(1, 6),rep(2, 15) )
)
europe_map <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
europe_map$iso_a2[europe_map$admin == "Greece"] <- "GR"
europe_map$iso_a2[europe_map$admin == "Cyprus"] <- "CY"
europe_map$iso_a2[europe_map$admin == "France"] <- "FR"
europe_map <- europe_map %>% filter(iso_a2 != "RU")
europe_clustered <- europe_map %>%
  left_join(cluster_data_harta, by = "iso_a2")
ggplot(data = europe_clustered) +
  geom_sf(aes(fill = as.factor(cluster)), color = "black") +
  scale_fill_manual(
    values = c( "1" = "#F8768D","2" = "#02BFC4"),
    na.value = "grey90",
    name = "Cluster"
  ) +
  coord_sf(
    xlim = c(-10, 35), 
    ylim = c(34, 71)     
  ) +
  theme_minimal() +
  labs(
    title = "Clasificarea țărilor europene în două clustere"
  )


scores <- pcr_data[, c("CP1", "CP2", "CP3")]
d_mah <- mahalanobis(scores, colMeans(scores), cov(scores))
pcr_data$mahalanobis <- d_mah
boxplot(pcr_data$mahalanobis, main = "Distanțe Mahalanobis")

#------------------------------------------------------------------------------------------------------------
#Regresie PCR
#-----------------------------------------------------------------------------------------------------------------
pca_scores_r <- data.frame(acp$scores[, 1:3])
colnames(pca_scores_r) <- c("CP1", "CP2", "CP3")
pcr_data <- data.frame(
  Y = data1$X5,
  CP1 = pca_scores_r$CP1,
  CP2 = pca_scores_r$CP2,
  CP3 = pca_scores_r$CP3
)
pcr_data <- pca_scores_r
pcr_data$Y <- data1$X5
pcr_data$cod <- rownames(date_bun)

lm_model1 <- lm(Y ~ CP1+ CP2+ CP3 , data = pcr_data) 
summary(lm_model1)
lm_model2 <- lm(Y ~ CP1+ CP2 , data = pcr_data) 
summary(lm_model2)
lm_model3 <- lm(Y ~ CP1+ CP3 , data = pcr_data) 
summary(lm_model3)

AIC(lm_model1, lm_model2, lm_model3)
BIC(lm_model1, lm_model2, lm_model3)

(plot(lm_model1))
(plot(lm_model2))
(plot(lm_model3))

pcr_data$Predicted <- predict(lm_model1)
pcr_data$Residuals <- resid(lm_model1)

pred_table <- data.frame(
  Observat = pcr_data$Y,
  Prezis = pcr_data$Predicted,
  Eroare = pcr_data$Y - pcr_data$Predicted
)
print(pred_table)

#graficul regresiei
ggplot(pcr_data, aes(x = Y, y = Predicted)) +
  geom_point(size = 3, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Valori Observate vs. Prezise",
    x = "Valori Observate (X5)",
    y = "Valori Prezise (model regresie)"
  ) +
  theme_minimal()

#histograma reziduurilor
ggplot(pcr_data, aes(x = Residuals)) +
  geom_histogram(bins = 10, fill = "skyblue", color = "black") +
  labs(title = "Distribuția Reziduurilor", x = "Reziduuri", y = "Frecvență") +
  theme_minimal()

#Testul Shapiro-Wilk
shapiro_test <- shapiro.test(pcr_data$Residuals)
print(shapiro_test)
durbinWatsonTest(lm_model)
bptest(lm_model)

#Grafic 3D 
plot_ly(
  data = pcr_data,
  x = ~CP1,
  y = ~CP2,
  z = ~CP3,
  text = rownames(date_bun), 
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers+text",
  marker = list(size = 5, color = ~Y, colorscale = "Blues"),
  textposition = "top middle"
) %>%
  layout(
    title = "Scoruri PCA colorate după Y (X11)",
    scene = list(
      xaxis = list(title = "CP1"),
      yaxis = list(title = "CP2"),
      zaxis = list(title = "CP3")
    )
  )
#Grafic 3D
grid <- expand.grid(
  CP1 = seq(min(pcr_data$CP1), max(pcr_data$CP1), length.out = 30),
  CP3 = seq(min(pcr_data$CP3), max(pcr_data$CP3), length.out = 30)
)
grid$CP2 <- mean(pcr_data$CP2)
grid$Y   <- predict(lm_model, newdata = grid)
z_matrix <- matrix(grid$Y, nrow = 30, ncol = 30)
plot_ly() %>%
  add_trace(
    data  = pcr_data,
    x     = ~CP1, y = ~CP3, z = ~Y,
    type  = "scatter3d", mode = "markers",
    marker = list(size = 5,
                  color = ~Y,
                  colorscale = "Blues",
                  showscale = FALSE),
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "CP1: %{x:.2f}<br>",
      "CP3: %{y:.2f}<br>",
      "Y  : %{z:.2f}<extra></extra>"
    ),
    text  = ~cod,
    name  = "Țări"
  ) %>%

  add_trace(
    data  = pcr_data,
    x     = ~CP1, y = ~CP3, z = ~Y,
    type  = "scatter3d", mode = "text",
    text  = ~cod,
    textfont = list(color = "black", size = 10),
    hoverinfo = "none",
    showlegend = FALSE
  ) %>%
  add_surface(
    x = unique(grid$CP1),
    y = unique(grid$CP3),
    z = z_matrix,
    opacity = 0.5,
    showscale = FALSE,
    name = "Plan regresie"
  ) %>%
  layout(
    title = "Model de regresie (lm) în spațiul PCA",
    scene = list(
      xaxis = list(title = "CP1"),
      yaxis = list(title = "CP3"),
      zaxis = list(title = "Y (X5)")
    )
  )

