"reduce2021" <- 
function(data, range)
{
	etendue <- diff(range)
	x <- 1/etendue
	y <- sweep(data, 2, FUN = "*", x)
	y <- sweep(y, 2, FUN = "+", (0 - range[1,  ]) * x)
	y
}

