## YUNUS IRMAK B1705.090008

##Get sequences

#<<<<<<<<<<<<<-FOR USER INPUT->>>>>>>>>>>>>>>>>>
#X <- readline(prompt="Enter the first sequence: ")
#Y <- readline(prompt="Enter the second sequence: ")
#<<<<<<<<<<<<<-FOR USER INPUT->>>>>>>>>>>>>>>>>>

X <- 'PAWHEAE'
Y <- 'HEAGAWGHEE'
seq.x <- unlist(strsplit(X, ''))
seq.y <- unlist(strsplit(Y, ''))

seq.x <- c(0,seq.x)
seq.y <- c(0,seq.y)

match <- 8
mismatch <- -8
indel <- -8

## creating the score matrix
xlength <- length(seq.x)
ylength <- length(seq.y)
score <- matrix(NA, xlength, ylength)
score[,1] <- sapply(1:xlength-1, function(x) x * 0)
score[1,] <- sapply(1:ylength-1, function(x) x * 0)
dimnames(score) <- list(seq.x[1:xlength], seq.y[1:ylength])

## The dynamic programming, global alignment recursion
for (i in 2:length(seq.x)) {
  for (j in 2:length(seq.y)){
    # seq.x[i] , seq.y[j] are aligned
    if ( seq.x[i] == seq.y[j]) {
      score[i,j] <- score[i-1, j-1] + match
    } else {
      score[i,j] <- score[i-1, j-1] + mismatch
      if (score[i,j]<0){
        score[i,j] <- 0
      }
    }
    # seq.x[i] aligned to -
    sc <- score[i-1,j] + indel
    if (score[i,j]<0){
      score[i,j] <- 0
    }
    if (sc > score[i,j])
      score[i,j] = sc
    # seq.y[j] aligned to -
    sc <- score[i,j-1] + indel
    if (score[i,j]<0){
      score[i,j] <- 0
    }
    if (sc > score[i,j])
      score[i,j] = sc
    if (score[i,j]<0){
      score[i,j] <- 0
    }
  }
}

## Traceback
i <- length(seq.x)
j <- length(seq.y)
ax <- character()
ay <- character()
while (i > 2 && j >2){
  ## case 1: best was seq.x[i] aligned to seq.y[j]
  sc <- score[i-1,j-1]
  if (seq.x[i] == seq.y[j]) {
    sc <- sc + match
  } else  {
    sc <- sc + mismatch
  }
  if (sc == score[i,j]) {
    ax <- c(seq.x[i], ax)
    ay <- c(seq.y[j], ay)
    i <- i -1
    j <- j-1
    next
  }
  ## case 2: best was seq.x[i] aligned to -
  if ((score[i-1,j] + indel) == score[i,j]) {
    ax <- c(seq.x[i], ax)
    ay <- c("-", ay)
    i <- i-1
    next
  }
  ## case 3: best was seq.y[j] aligned to -
  if ((score[i,j-1] + indel) == score[i,j]) {
    ax <- c("-", ax)
    ay <- c(seq.y[j], ay)
    j <- j-1
    next
  }
}

cat ("First Sequence: ", X,"\n")
cat ("Second Sequence: ", Y,"\n")
cat ("Scoring System: ", match, " for match; ", mismatch, " for mismatch; ", indel, " for gap", "\n\n")

cat ("Dynamic Programming Matrix:\n")
print (score)

cat ("\nAlignment:\n")
cat (paste(ax, collapse=''), "\n")
cat (paste(ay, collapse=''),"\n\n")
cat ("Optimum alignment score: ", score[length(score)],"\n")