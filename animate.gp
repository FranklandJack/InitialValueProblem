# Set a title.
set title "Cahn-Hilliard"

# Remove border.
unset border

# Remove x and y ticks.
unset xtics
unset ytics

# Remove color gradient.
#unset colorbox

# Set spin colours.
#set palette defined (0 "red", 1 "green", 2 "blue")

# This stops an error when all spins are the same.
set cbrange[-1:1]

# Do the plot.
p filename matrix with image
while(1){
	replot
	pause 0.3
}