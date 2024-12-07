# Parameters
nident <- 600                # Initial homology score
prot_genomic_len <- 1000     # Expected genomic length covered by the protein
maxIntronSize <- 500         # Maximum expected intron size

# Generate a sequence of class_length values
class_length <- seq(500, 3000, by = 1)

# Calculate Adjusted_score using the continuous penalty function
#Penalty <- 0.05 * pmax(0, (class_length - prot_genomic_len) / maxIntronSize)
class_lengthOK=pmax(prot_genomic_len, class_length)
pc_penalty=0.01

Penalty<-pmax(-1, pc_penalty*(1-pmax(1,exp((class_length - prot_genomic_len)/maxIntronSize))))
Adjusted_score <- nident *(1+Penalty)
Penalty
print (Bonus[1])
print (Penalty[1])
# Plot the scoring function
plot(class_length, Adjusted_score, type = "l", col = "blue", lwd = 2,
     xlab = "Class Length",
     ylab = "Adjusted Score",
     main = "Adjusted Score vs. Class Length (Continuous Penalty)")
# Add a horizontal line at nident
abline(h = nident, col = "red", lty = 2)

# Add vertical line at prot_genomic_len
abline(v = prot_genomic_len, col = "green", lty = 2)

# Add legend
legend("topright", legend = c("Adjusted Score", "nident", "prot_genomic_len"),
       col = c("blue", "red", "green"), lty = c(1, 2, 2), lwd = c(2, 1, 1))
log10(10)
