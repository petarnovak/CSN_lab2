#1: Introduction
write_table <- function(language,file,summarytab.df) {
   degree_sequence = read.table(file, header = FALSE)
   summarytab.df <- rbind(summarytab.df,data.frame("Language" = language, 
                                                   "N" = length(degree_sequence$V1), 
                                                   "Maximum degree" = max(degree_sequence$V1),
                                                   "M/N"=sum(degree_sequence$V1)/length(degree_sequence$V1),
                                                   "N/M"=length(degree_sequence$V1)/sum(degree_sequence$V1)))
   return(summarytab.df)
   }

source = read.table("list_in.txt", 
         header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
         as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
        )

summarytab.df <- data.frame("Language"=character(), "N"=numeric(),"Maximum degree"=numeric(), 
                            "M/N"=numeric(),"N/M"=numeric())

for (x in 1:nrow(source)) {
  summarytab.df <- write_table(source$language[x], source$file[x],summarytab.df)
}
summarytab.df
