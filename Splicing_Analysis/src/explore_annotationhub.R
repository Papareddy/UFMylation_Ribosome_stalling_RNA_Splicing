
library(AnnotationHub)

print("Connecting to AnnotationHub...")
ah <- AnnotationHub()

# Query for CisBP RNA
print("Querying for 'CisBP' and 'RNA'...")
res <- query(ah, c("CisBP", "RNA", "Homo sapiens"))
print(res)

# Ray2013?
print("Querying for 'Ray2013'...")
res2 <- query(ah, c("Ray2013"))
print(res2)

# If found, how to download?
# Usually: m <- res[[id]]
# Then check format.

# ATtRACT?
print("Querying for 'ATtRACT'...")
res3 <- query(ah, c("ATtRACT"))
print(res3)
