
################################################################
### Program to translate the 09 curves to the standard format:
### cal. BP   rc BP    std. error,
### in cal. BP increasing order, with no commas.
#################################################################


### Read the terrestrial and Marine curves files as distributed in
### http://www.radiocarbon.org/IntCal09.htm :

cc.terr <- read.csv( "intcal09.14c", header=FALSE, skip=11)
cc.marine <- read.csv( "marine09.14c", header=FALSE, skip=11)

### Write them in the desired format:

write.table( cc.terr[nrow(cc.terr):1,1:3], file="3Col_intcal09.14C", row.names=FALSE, col.names=FALSE)

write.table( cc.marine[nrow(cc.marine):1,1:3], file="3Col_marine09.14C", row.names=FALSE, col.names=FALSE)

