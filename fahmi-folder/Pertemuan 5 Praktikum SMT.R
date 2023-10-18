
library(readxl)
dataku <- read_excel("D:/ASPRAK/SMT/PERTEMUAN 5/Data Pertemuan 5.xlsx")

# Menghitung jumlah missing value pada setiap variabel
missing_count <- sapply(dataku, function(x) sum(is.na(x)))

# Menghitung persentase missing value pada setiap variabel
total_rows <- nrow(dataku)
missing_percentage <- (missing_count / total_rows) * 100

# Menampilkan hasil
result <- data.frame(Var = names(missing_count), MissingCount = missing_count, MissingPercentage = missing_percentage)
result

# Mengisi missing value dengan rata-rata
dataku <- dataku %>%
  mutate_if(is.numeric, ~ifelse(is.na(.), mean(., na.rm = TRUE), .))
sapply(dataku, function(x) sum(is.na(x)))

#####
# Menguji adanya outlier
#a. standardisasi
data_outlier <- read_excel("D:/ASPRAK/SMT/PERTEMUAN 5/Data Pertemuan 5.xlsx", 
                               sheet = "Outlier")
standardized_data <- data_outlier %>%
  mutate(
    zusia = scale(usia),
    zberat = scale(berat),
    ztinggi = scale(tinggi),
    zincome = scale(income),
    zjamkerja = scale(jamker),
    zolahraga = scale(olhraga)
  )
standardized_data

#b. Scatter plot
# Membuat scatterplot
library(ggplot2)
ggplot(data_outlier, aes(x = income, y = usia)) +
  geom_point() +
  labs(x = "Income", y = "Usia") +
  ggtitle("Scatterplot Usia vs. Income")

#c. boxplot
ggplot(data_outlier, aes(y = usia)) +
  geom_boxplot(fill = "skyblue", color = "blue") +  # Atur warna kotak dan garis tepi
  labs(y = "Usia") +
  ggtitle("Boxplot Usia") +
  theme_minimal()
