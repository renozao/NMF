# NMF package
#
# Helper code to allow mixing grid/base graphics.
# This code is only loaded with the explicit request of the user, 
# either via option 'grid.patch' or environment variable R_PACKAGE_NMF_GRID_PATCH.  
#
# The functions in this file were adapted from the grid package, which is 
# under the following GPL license:
#
# R : A Computer Language for Statistical Data Analysis
# Copyright (C) 2001-3 Paul Murrell
#               2003 The R Core Team
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, a copy is available at
# http://www.r-project.org/Licenses/
# 
# Author: Renaud Gaujoux
# Created: Sep 16, 2013
###############################################################################

# This is essentially where the patch lies: not calling L_gridDirty 
grid.Call <- function (fnname, ...) 
{
	#.Call(L_gridDirty)
	.Call(fnname, ..., PACKAGE = "grid")
}

# One has to test for nullity since not using L_gridDirty means potentially 
# returning a NULL viewport 
current.viewport <- function()
{
    cv <- grid.Call(grid:::L_currentViewport)
    if( !is.null(cv) ) grid:::vpFromPushedvp(cv)
}

# same thing here: call patched current.viewport and 
# check for nullity
current.vpPath <- function(){
 	names <- NULL
 	pvp <- current.viewport()
 	if( is.null(pvp) ) return(NULL)
 	while ( !is.null(pvp) && !grid:::rootVP(pvp)) {
 		names <- c(names, pvp$name)
 		pvp <- pvp$parent
 	}
 	if (!is.null(names)) 
 		grid:::vpPathFromVector(rev(names))
 	else names	
}
