#1: Introduction
write_table <- function(language,file,lang.df) {
   degree_sequence = read.table(file, header = FALSE)
   lang.df <- rbind(lang.df,data.frame(language, length(degree_sequence$V1), max(degree_sequence$V1),
                                       sum(degree_sequence$V1)/length(degree_sequence$V1),
                                       length(degree_sequence$V1)/sum(degree_sequence$V1)))
   return(lang.df)
   }

source = read.table("list_in.txt", 
         header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
         as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
        )

lang.df <- data.frame()

for (x in 1:nrow(source)) {
  lang.df <- write_table(source$language[x], source$file[x],lang.df)
}
colnames(lang.df) <- c("Language", "N", "Maximum degree", "M/N", "N/M")
#Table with all languages and some important values
lang.df

#2: Visualization
# Function to plot the degree sequence of a given language
degree_seq <- function(language){
  path = paste("./data/",language,"_in-degree_sequence.txt",sep="")
  degree_sequence = read.table(path, header = FALSE)
}

#Arabic

#English
degree_sequence = read.table("./data/English_in-degree_sequence.txt", header = FALSE)
degree_spectrum = table(degree_sequence)
degree_sequence
plot(rownames(degree_sequence),degree_sequence[,1])
barplot(degree_spectrum, main = "English",
        xlab = "degree", ylab = "number of vertices")
barplot(degree_spectrum, main = "English",
        xlab = "degree", ylab = "number of vertices", log = "xy")
