"square" <- 
function(rayon = 1, centre.x = 0, centre.y = 0)
{
	x <- centre.x + c( - rayon,  - rayon, rayon, rayon)
	y <- centre.y + c(rayon,  - rayon,  - rayon, rayon)
	list(x = x, y = y)
}

