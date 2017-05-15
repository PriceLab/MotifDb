# In order to correctly format my matrix for the parsePWMs function, had to do the following:
text <- gsub("[A,C,G,T]\\s+\\[\\w?","",text)
text <- gsub("\\s+\\]","",text)
text <- gsub("\\s+","\t",text)

# This is because the file had the following issues:
## Row names were there
## Brackets were there
## Some opening brackets had spaces and others didn't
## Entries were space-delimited

## Best now to change to input file OR change the processing steps?


# Also, titles have duplicates in the metadata table...which columns shouldn't be duplicated?
# id
# ma.name

# writing data to file:
write.table(tbl.md, file="md-raw.tsv", sep = "\t", row.names = FALSE)


