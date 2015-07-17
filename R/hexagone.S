"hexagone" <- 
function(rayon = 1, offset.angle = 0, centre.x = 0, centre.y = 0)
{
	angles <- seq(0, 2 * pi, pi/3) + offset.angle
	x <- cos(angles) * rayon + centre.x
	y <- sin(angles) * rayon + centre.y
	res <- list(x = x, y = y)
	res
}

