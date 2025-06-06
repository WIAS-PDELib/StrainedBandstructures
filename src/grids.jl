function bimetal_strip3D(; material_border = 0.5, scale = [1,1,1], anisotropy = [1,1,2], reflevel = 1, maxvol = prod(scale./anisotropy)/8)

    @info "Generating bimetal grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=TetGen)
    p1=point!(builder,0,0,0)
    p2=point!(builder,scale[1],0,0)
    p3=point!(builder,scale[1],scale[2],0)
    p4=point!(builder,0,scale[2],0)
    p5=point!(builder,0,0,scale[3])
    p6=point!(builder,scale[1],0,scale[3])
    p7=point!(builder,scale[1],scale[2],scale[3])
    p8=point!(builder,0,scale[2],scale[3])

    p9=point!(builder,0,material_border*scale[2],0)
    p10=point!(builder,scale[1],material_border*scale[2],0)
    p11=point!(builder,scale[1],material_border*scale[2],scale[3])
    p12=point!(builder,0,material_border*scale[2],scale[3])

    facetregion!(builder,1) # bottom (material A) = core front
    facet!(builder,p1 ,p2 ,p10 ,p9)
    facet!(builder,p9 ,p10 ,p3 ,p4)

    facetregion!(builder,2) # shell boundary
    facet!(builder,p6 ,p2 ,p1 ,p5)
    facet!(builder,p11 ,p6 ,p5 ,p12)
    facet!(builder,p2 ,p6 ,p11 ,p10)
    facet!(builder,p1 ,p9 ,p12 ,p5)

    facetregion!(builder,3) # stressor boundary
    facet!(builder,p7 ,p11 ,p12 ,p8)
    facet!(builder,p10 ,p11 ,p7 ,p3)
    facet!(builder,p3 ,p7 ,p8 ,p4)
    facet!(builder,p9 ,p4 ,p8 ,p12)

    facetregion!(builder,4) # interior facet to split regions
    facet!(builder,p9 ,p10 ,p11 ,p12)

    cellregion!(builder,1)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[1],0.5*material_border*scale[2],0.5*scale[3]))

    cellregion!(builder,2)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[1],0.5*(material_border+1)*scale[2],0.5*scale[3]))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 3, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)
    return xgrid

end


function bimetal_strip3D_middle_layer(; material_border = 0.5, scale = [1,1,1], anisotropy = [1,1,2], reflevel = 1, maxvol1 = prod(scale./anisotropy)/8, maxvol2 = prod(scale./anisotropy))

    @info "Generating bimetal sandwich grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=TetGen)
    y = 5
    p1=point!(builder,0,0,0)
    p2=point!(builder,scale[1],0,0)
    p3=point!(builder,scale[1],material_border*scale[2]-y,0)
    p4=point!(builder,scale[1],material_border*scale[2],0)
    p5=point!(builder,scale[1],material_border*scale[2]+y,0)
    p6=point!(builder,scale[1],scale[2],0)
    p7=point!(builder,0,scale[2],0)
    p8=point!(builder,0,material_border*scale[2]+y,0)
    p9=point!(builder,0,material_border*scale[2],0)
    p10=point!(builder,0,material_border*scale[2]-y,0)

    p11=point!(builder,0,0,scale[3])
    p12=point!(builder,scale[1],0,scale[3])
    p13=point!(builder,scale[1],material_border*scale[2]-y,scale[3])
    p14=point!(builder,scale[1],material_border*scale[2],scale[3])
    p15=point!(builder,scale[1],material_border*scale[2]+y,scale[3])
    p16=point!(builder,scale[1],scale[2],scale[3])
    p17=point!(builder,0,scale[2],scale[3])
    p18=point!(builder,0,material_border*scale[2]+y,scale[3])
    p19=point!(builder,0,material_border*scale[2],scale[3])
    p20=point!(builder,0,material_border*scale[2]-y,scale[3])

    facetregion!(builder,1) # core boundary (bottom, top, sides)
    facet!(builder,p1,p2,p3,p10)
    facetregion!(builder,12)
    facet!(builder,p11,p12,p13,p20)
    facetregion!(builder,13)
    facet!(builder,p1,p2,p12,p11)
    facet!(builder,p2,p3,p13,p12)
    facet!(builder,p1,p10,p20,p11)

    facetregion!(builder,2) # middle layer 1 (bottom, top, sides)
    facet!(builder,p3,p4,p9,p10)
    facet!(builder,p13,p14,p19,p20)
    facet!(builder,p3,p10,p20,p13)
    facet!(builder,p3,p4,p14,p13)
    facet!(builder,p9,p10,p20,p19)

    facetregion!(builder,3) # middle layer 2 (bottom, top, sides)
    facet!(builder,p4,p5,p8,p9)
    facet!(builder,p14,p15,p18,p19)
    facet!(builder,p5,p8,p18,p15)
    facet!(builder,p4,p5,p15,p14)
    facet!(builder,p8,p9,p19,p18)

    facetregion!(builder,4) # stressor boundary (bottom, top, sides)
    facet!(builder,p5,p6,p7,p8)
    facet!(builder,p15,p16,p17,p18)
    facet!(builder,p6,p7,p17,p16)
    facet!(builder,p5,p6,p16,p15)
    facet!(builder,p7,p8,p18,p17)

    facetregion!(builder,5) # interior facet to split regions
    facet!(builder,p4,p9,p19,p14)

    cellregion!(builder,1)
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*scale[1],0.5*(material_border*scale[2]-y),0.5*scale[3]))
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*scale[1],material_border*scale[2]-y/2,0.5*scale[3]))

    cellregion!(builder,2)
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*scale[1],0.5*((material_border+1)*scale[2]+y),0.5*scale[3]))
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*scale[1],material_border*scale[2]+y/2,0.5*scale[3]))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 3, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)
    return xgrid

end


function bimetal_strip2D(; material_border = 0.5, scale = [1,1], anisotropy = [1,1], reflevel = 1, maxvol = prod(scale./anisotropy)/4)

    @info "Generating 2d bimetal grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=Triangulate)
    p1=point!(builder,0,0)
    p2=point!(builder,scale[2],0)
    p3=point!(builder,scale[2],scale[1])
    p4=point!(builder,0,scale[1])

    p5=point!(builder,0,material_border*scale[1])
    p6=point!(builder,scale[2],material_border*scale[1])


    facetregion!(builder,1) # left (material A)
    facet!(builder,p5 ,p1)
    facetregion!(builder,11) # left (material B)
    facet!(builder,p4,p5)

    facetregion!(builder,2) # bottom
    facet!(builder,p1 ,p2)

    facetregion!(builder,3) # top
    facet!(builder,p3, p4)

    facetregion!(builder,4) # right
    facet!(builder,p2, p6)
    facet!(builder,p6, p3)

    facetregion!(builder,99) # interior facet to split regions
    facet!(builder,p5, p6)

    cellregion!(builder,1)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[2],0.5*material_border*scale[1]))

    cellregion!(builder,2)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[2],0.5*(material_border+1)*scale[1]))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 2, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)
    return xgrid

end


function condensator3D(; scale=[50,50,50], d = 5, nrefs = 1, maxvol1 = prod(scale)/4, maxvol2 = scale[1]*scale[2]*d/4)

    builder=SimplexGridBuilder(Generator=TetGen)

    @info "Generating 3D condensator grid for scale = $scale and middle layer width = $d."

	X = scale[1]
	Y = scale[2]
	Z = scale[3]

    # side at y = 0
    p1=point!(builder,0,0,0)
    p2=point!(builder,X,0,0)
    p3=point!(builder,X,0,Z)
    p4=point!(builder,X,0,Z+d)
    p5=point!(builder,X,0,2*Z+d)
    p6=point!(builder,0,0,2*Z+d)
    p7=point!(builder,0,0,Z+d)
    p8=point!(builder,0,0,Z)

    # side at y = Y
    p9=point!(builder,0,Y,0)
    p10=point!(builder,X,Y,0)
    p11=point!(builder,X,Y,Z)
    p12=point!(builder,X,Y,Z+d)
    p13=point!(builder,X,Y,2*Z+d)
    p14=point!(builder,0,Y,2*Z+d)
    p15=point!(builder,0,Y,Z+d)
    p16=point!(builder,0,Y,Z)

    facetregion!(builder,1) # bottom of bottom cube
    facet!(builder,p1,p9,p10,p2)
    facetregion!(builder,2) # front, back, left and right boundary of bottom cube
    facet!(builder,p1,p2,p3,p8)
    facet!(builder,p10,p9,p16,p11)
    facet!(builder,p2,p10,p11,p3)
    facet!(builder,p9,p1,p8,p16)

    facetregion!(builder,3) # gap boundary front, back
    facet!(builder,p8,p3,p4,p7)
    facet!(builder,p11,p16,p15,p12)
    facetregion!(builder,4) # gap boundary left and right
    facet!(builder,p3,p11,p12,p4)
    facet!(builder,p16,p8,p7,p15)

    facetregion!(builder,5) # boundary of top cube (front, back, left, right)
    facet!(builder,p7,p4,p5,p6)
    facet!(builder,p12,p15,p14,p13)
    facet!(builder,p4,p12,p13,p5)
    facet!(builder,p15,p7,p6,p14)
    facetregion!(builder,6) # boundary of top cube (top)
    facet!(builder,p6,p5,p13,p14)

    #facetregion!(builder,7) # interior facets to split materials
    facet!(builder,p8,p3,p11,p16)
    facet!(builder,p7,p4,p12,p15)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*X,0.5*Y,0.5*Z))
    regionpoint!(builder,(0.5*X,0.5*Y,1.5*Z+d))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*X,0.5*Y,Z+0.5*d))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,nrefs)

    return xgrid
end

function condensator3D_tensorgrid(; scale = [50,50,50], d = scale[3]/5, nrefs = 2)

    @info "Generating 3D condensator grid for a cuboid with dimensions ($(scale[1]),$(scale[2]),$(2*scale[1]+d)) and middle layer of width $d."

    X = scale[1]
    Y = scale[2]
    Z = 2*scale[3] + d
    npts = 2*(nrefs+1)

    XX = LinRange(0,scale[1],npts)
    YY = LinRange(0,scale[2],npts)

    xgrid = simplexgrid(XX,YY)
    xgrid_cross_section = deepcopy(xgrid)

    z_levels_nonuniform = Array{Float64,1}(LinRange(0,scale[3]-d,npts))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3]-d,scale[3],1+Int(npts/2))[3:end]))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3],scale[3]+d,1+Int(npts/2))[2:end-1]))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3]+d,scale[3]+2*d,1+Int(npts/2))[1:end-2]))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3]+2*d,Z,npts)))

    xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=4, top_offset=5)
    xgrid = uniform_refine(xgrid,nrefs-1)
    # the offsets lead to the following boundary regions:
    # 1 - 4  = side core
    # 5      = bottom core
    # 6      = tope core
    # 7 - 10 = side stressor

    cellmask!(xgrid,[0,0,scale[3]],[X,Y,scale[3]+d],2)
	bfacemask!(xgrid,[0,0,scale[3]],[X,0,scale[3]+d],7)
	bfacemask!(xgrid,[0,Y,scale[3]],[X,Y,scale[3]+d],8)
    bfacemask!(xgrid,[0,0,scale[3]],[0,Y,scale[3]+d],9)
	bfacemask!(xgrid,[X,0,scale[3]],[X,Y,scale[3]+d],10)

    return xgrid, xgrid_cross_section
end

function condensator3D_tensorgrid!(; scale=[50,50,50], d=10, nrefs=0, dx=0.5,
	stressor_cell_per=10)

	X = scale[1]
    Y = scale[2]
    Z0 = scale[3]
	Z = 2*scale[3] + d

    @info "Generating 3D condensator grid for a cuboid with dimensions
			($(X),$(Y),$(Z)) and middle layer of width $d."

    XX = LinRange(0,X,Int(round(X/dx))+1)
    YY = LinRange(0,Y,Int(round(Y/dx))+1)

	xygrid = simplexgrid(XX,YY)
    grid_cross_section = deepcopy(xygrid)

	z_levels_uniform = LinRange(0,Z0,Int(round(Z0/dx))+1)
	bottom_grid = simplexgrid(xygrid, z_levels_uniform, bot_offset=4, top_offset=5)

	z_levels_uniform = LinRange(Z0,Z0+d,Int(round(d/dx))+1)
	middle_grid = simplexgrid(xygrid, z_levels_uniform, bot_offset=4, top_offset=5)

	z_levels_uniform = LinRange(Z0+d,Z,Int(round(Z0/dx))+1)
	top_grid = simplexgrid(xygrid, z_levels_uniform, bot_offset=4, top_offset=5)

	cell_num = Int(round(stressor_cell_per/100*length(middle_grid[CellRegions])/6))
	cell_indices = shuffle(1:6:length(middle_grid[CellRegions]))[1:cell_num]
	for (~,j) in enumerate(cell_indices)
		middle_grid[CellRegions][j:j+5] .= 2
	end

	grid = glue(bottom_grid,middle_grid)
	grid = glue(grid,top_grid)

    # the offsets lead to the following boundary regions (anti-clockwise indexing):
    # 1 - 4  = side core
    # 5      = bottom core
    # 6      = tope core
    # 7 - 10 = side stressor
	bfacemask!(grid,[0,0,Z0],[X,0,Z0+d],7)
	bfacemask!(grid,[X,0,Z0],[X,Y,Z0+d],8)
	bfacemask!(grid,[0,Y,Z0],[X,Y,Z0+d],9)
	bfacemask!(grid,[0,0,Z0],[0,Y,Z0+d],10)

    return grid, grid_cross_section
end

function condensator2D(; A = 50, B = 100, d = 5, reflevel = 1, maxvol1 = B*d/4, maxvol2 = B*d/4)

    builder=SimplexGridBuilder(Generator=Triangulate)

    @info "Generating 2d condensator grid for A = $A, B = $B, d = $d"

    p1=point!(builder,0,0)
    p2=point!(builder,B,0)
    p3=point!(builder,B,A)
    p4=point!(builder,B,A+d)
    p5=point!(builder,B,2*A+d)
    p6=point!(builder,0,2*A+d)
    p7=point!(builder,0,A+d)
    p8=point!(builder,0,A)

    facetregion!(builder,1)
    facet!(builder,p1,p2) # bottom of bottom plate

    facetregion!(builder,2)
    facet!(builder,p2,p3) # right boundary of bottom plate
    facet!(builder,p3,p4) # right gap boundary
    facet!(builder,p4,p5) # right boundary of top plate

    facetregion!(builder,3)
    facet!(builder,p5,p6) # top boundary of top plate

    facetregion!(builder,4)
    facet!(builder,p6,p7) # left boundary of top plate
    facet!(builder,p7,p8) # left gap boundary
    facet!(builder,p8,p1) # left boundary of bottom plate

    # interior facets to split materials
	facetregion!(builder,6)
    facet!(builder,p3,p8)
    facet!(builder,p4,p7)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*B,0.5*A))
    regionpoint!(builder,(0.5*B,1.5*A+d))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*B,A+0.5*d))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function condensator2D_tensorgrid(; scale = [50,50], d = 10, nrefs = 2)

    @info "Generating 2D condensator grid for a rectangle with dimensions ($(scale[1]),$(2*scale[2]+d)) and middle layer of width $d."

    X = scale[1]
    Y = 2*scale[2] + d
    npts = 2*(nrefs+1)

    XX = LinRange(0,scale[1],npts)
    YY = Array{Float64,1}(LinRange(0,scale[2]-d,npts))
    append!(YY, Array{Float64,1}(LinRange(scale[2]-d,scale[2],1+Int(npts/2))[3:end]))
    append!(YY, Array{Float64,1}(LinRange(scale[2],scale[2]+d,1+Int(npts/2))[2:end-1]))
    append!(YY, Array{Float64,1}(LinRange(scale[2]+d,scale[2]+2*d,1+Int(npts/2))[1:end-2]))
    append!(YY, Array{Float64,1}(LinRange(scale[2]+2*d,Y,npts)))

	xgrid = simplexgrid(XX,YY)
    xgrid = uniform_refine(xgrid,nrefs-1)
    # the offsets lead to the following boundary regions:
    # 1 = bottom bulk
    # 2 = right bulk
    # 3 = top bulk
    # 4 = left bulk
	# 5 = right side stressor
	# 6 = left side stressor

    cellmask!(xgrid,[0,scale[2]],[X,scale[2]+d],2)
    bfacemask!(xgrid,[X,scale[2]],[X,scale[2]+d],5)
    bfacemask!(xgrid,[0,scale[2]],[0,scale[2]+d],6)

    return xgrid
end


function nonpolarquantumwell3D(; X = 80, Y = 50, Z = 20, V = 8, d = 2, w = 3, reflevel = 1, maxvol1 = 1/4*(Y*Z*(X-V)/2), maxvol2 = 1/4*(Y*Z*V))

    builder=SimplexGridBuilder(Generator=TetGen)

    @info "Generating 3d quantum well grid for X = $X, Y = $Y, Z = $Z, V = $V, d = $d, w = $w"

    # bottom side at Y = 0
    p1=point!(builder,0,0,0)
    p2=point!(builder,0,0,Z)
    p3=point!(builder,(X-V)/2,0,Z)
    p4=point!(builder,(X-V)/2+V,0,Z)
    p5=point!(builder,X,0,Z)
    p6=point!(builder,X,0,0)
    p7=point!(builder,(X-V)/2+V,0,0)
    p8=point!(builder,(X-V)/2,0,0)
    p9=point!(builder,(X-V)/2+V,0,(Z+w)/2)
    p10=point!(builder,(X-V)/2+V-d,0,(Z+w)/2)
    p11=point!(builder,(X-V)/2+V-d,0,(Z-w)/2)
    p12=point!(builder,(X-V)/2+V,0,(Z-w)/2)

    # top side at Y
    p13=point!(builder,0,Y,0)
    p14=point!(builder,0,Y,Z)
    p15=point!(builder,(X-V)/2,Y,Z)
    p16=point!(builder,(X-V)/2+V,Y,Z)
    p17=point!(builder,X,Y,Z)
    p18=point!(builder,X,Y,0)
    p19=point!(builder,(X-V)/2+V,Y,0)
    p20=point!(builder,(X-V)/2,Y,0)
    p21=point!(builder,(X-V)/2+V,Y,(Z+w)/2)
    p22=point!(builder,(X-V)/2+V-d,Y,(Z+w)/2)
    p23=point!(builder,(X-V)/2+V-d,Y,(Z-w)/2)
    p24=point!(builder,(X-V)/2+V,Y,(Z-w)/2)

    facetregion!(builder,1) # sides of region 1 cube (bottom, top, front, left, right)
    facet!(builder,p1,p2,p3,p8)
    facet!(builder,p13,p14,p15,p20)
    facet!(builder,p1,p2,p14,p13)
    facet!(builder,p1,p8,p20,p13)
    facet!(builder,p2,p3,p15,p14)

    facetregion!(builder,2) # gap region
    # bottom
    facet!(builder,p3,p4,p9,p10)
    facet!(builder,p3,p10,p11,p8)
    facet!(builder,p8,p11,p12,p7)
    # top
    facet!(builder,p15,p16,p21,p22)
    facet!(builder,p15,p22,p23,p20)
    facet!(builder,p20,p23,p24,p19)
    # front
    facet!(builder,p3,p8,p20,p15)
    # back
    facet!(builder,p4,p9,p21,p16)
    facet!(builder,p9,p10,p22,p21)
    facet!(builder,p10,p11,p23,p22)
    facet!(builder,p11,p12,p24,p23)
    facet!(builder,p12,p7,p19,p24)
    # left
    facet!(builder,p7,p8,p20,p19)
    # right
    facet!(builder,p3,p4,p16,p15)

    facetregion!(builder,3) # sides of region 2 cube
    # bottom
    facet!(builder,p4,p5,p9)
    facet!(builder,p5,p9,p12,p6)
    facet!(builder,p6,p12,p7)
    #facet!(builder,p4,p5,p6,p9)
    #facet!(builder,p6,p9,p12,p7)
    facet!(builder,p9,p10,p11,p12)
    # top
    facet!(builder,p16,p17,p18,p21)
    facet!(builder,p18,p21,p24,p19)
    facet!(builder,p21,p22,p23,p24)
    # back
    facet!(builder,p5,p6,p18,p17)
    # left
    facet!(builder,p6,p7,p19,p18)
    # right
    facet!(builder,p4,p5,p17,p16)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*(X-V)/2,0.5*Y,0.5*Z))
    regionpoint!(builder,(X-0.5*(X-V)/2,0.5*Y,0.5*Z))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,((X-V)/2+0.5*d,0.5*Y,0.5*Z))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function nonpolarquantumwell2D(; X = 80, Z = 20, V = 8, d = 2, w = 3, reflevel = 1, maxvol1 = 1/4*(Z*(X-V)/2), maxvol2 = 1/4*(Z*V))

    builder=SimplexGridBuilder(Generator=Triangulate)

    @info "Generating 2d quantum well grid for X = $X, Z = $Z, V = $V, d = $d, w = $w"

    p1=point!(builder,0,0)
    p2=point!(builder,0,Z)
    p3=point!(builder,(X-V)/2,Z)
    p4=point!(builder,(X-V)/2+V,Z)
    p5=point!(builder,X,Z)
    p6=point!(builder,X,0)
    p7=point!(builder,(X-V)/2+V,0)
    p8=point!(builder,(X-V)/2,0)
    p9=point!(builder,(X-V)/2+V,(Z+w)/2)
    p10=point!(builder,(X-V)/2+V-d,(Z+w)/2)
    p11=point!(builder,(X-V)/2+V-d,(Z-w)/2)
    p12=point!(builder,(X-V)/2+V,(Z-w)/2)

    facetregion!(builder,1) # boundary of bottom plate (bottom, right, left)
    facet!(builder,p1,p2)
    facet!(builder,p2,p3)
    facet!(builder,p8,p1)

    facetregion!(builder,2) # gap boundary (right, left)
    facet!(builder,p3,p4)
    facet!(builder,p7,p8)

    facetregion!(builder,3) # interior facets to split materials
    facet!(builder,p3,p8)
    facet!(builder,p4,p9)
    facet!(builder,p9,p10)
    facet!(builder,p10,p11)
    facet!(builder,p11,p12)
    facet!(builder,p12,p7)

    facetregion!(builder,4) # boundary of top plate (left, top, right)
    facet!(builder,p4,p5)
    facet!(builder,p5,p6)
    facet!(builder,p6,p7)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*(X-V)/2),0.5*Z)
    regionpoint!(builder,(X-0.5*(X-V)/2),0.5*Z)

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,((X-V)/2+0.5*d),0.5*Z)

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function nanowire_grid(; scale = [1,1,1,1], anisotropy = [1,1,1,1], reflevel = 1, maxvol = prod(scale./anisotropy)/8)

    @info "Generating nanowire grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=TetGen)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ = scale[3]

    # bottom side at Z = 0
    p0=point!(builder,0,0,0)
    p1=point!(builder,d1,0,0)
    p2=point!(builder,d1/2,sqrt(3)/2*d1,0)
    p3=point!(builder,-d1/2,sqrt(3)/2*d1,0)
    p4=point!(builder,-d1,0,0)
    p5=point!(builder,-d1/2,-sqrt(3)/2*d1,0)
    p6=point!(builder,d1/2,-sqrt(3)/2*d1,0)
    p7=point!(builder,d2,0,0)
    p8=point!(builder,d2/2,sqrt(3)/2*d2,0)
    p9=point!(builder,-d2/2,sqrt(3)/2*d2,0)
    p10=point!(builder,-d2,0,0)
    p11=point!(builder,-d2/2,-sqrt(3)/2*d2,0)
    p12=point!(builder,d2/2,-sqrt(3)/2*d2,0)
    p13=point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ,0)
    p14=point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ,0)

    # top side at Z = scale[4]
    p15=point!(builder,d1,0,scale[4])
    p16=point!(builder,d1/2,sqrt(3)/2*d1,scale[4])
    p17=point!(builder,-d1/2,sqrt(3)/2*d1,scale[4])
    p18=point!(builder,-d1,0,scale[4])
    p19=point!(builder,-d1/2,-sqrt(3)/2*d1,scale[4])
    p20=point!(builder,d1/2,-sqrt(3)/2*d1,scale[4])
    p21=point!(builder,d2,0,scale[4])
    p22=point!(builder,d2/2,sqrt(3)/2*d2,scale[4])
    p23=point!(builder,-d2/2,sqrt(3)/2*d2,scale[4])
    p24=point!(builder,-d2,0,scale[4])
    p25=point!(builder,-d2/2,-sqrt(3)/2*d2,scale[4])
    p26=point!(builder,d2/2,-sqrt(3)/2*d2,scale[4])
    p27=point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ,scale[4])
    p28=point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ,scale[4])
    p29=point!(builder,0,0,scale[4])

    facetregion!(builder,1) # core region (bottom, sides, top)
    # bottom
    facet!(builder,p0,p1,p2)
    facet!(builder,p0,p2,p3)
    facet!(builder,p0,p3,p4)
    facet!(builder,p0,p4,p5)
    facet!(builder,p0,p5,p6)
    facet!(builder,p0,p6,p1)
    #facet!(builder,p1,p2,p3,p4)
    #facet!(builder,p4,p5,p6,p1)
    # sides
    facetregion!(builder,11) # core region (bottom, sides, top)
    facet!(builder,p1,p2,p16,p15)
    facet!(builder,p2,p3,p17,p16)
    facet!(builder,p3,p4,p18,p17)
    facet!(builder,p4,p5,p19,p18)
    facet!(builder,p5,p6,p20,p19)
    facet!(builder,p6,p1,p15,p20)
    # top
    facetregion!(builder,12) # core region (bottom, sides, top)
    facet!(builder,p29,p15,p16)
    facet!(builder,p29,p16,p17)
    facet!(builder,p29,p17,p18)
    facet!(builder,p29,p18,p19)
    facet!(builder,p29,p19,p20)
    facet!(builder,p29,p20,p15)
    #facet!(builder,p15,p16,p17,p18)
    #facet!(builder,p18,p19,p20,p15)

    facetregion!(builder,2) # shell region (bottom, sides, top)
    # bottom
    facet!(builder,p1,p7,p8,p2)
    facet!(builder,p2,p8,p9,p3)
    facet!(builder,p3,p9,p10,p4)
    facet!(builder,p4,p10,p11,p5)
    facet!(builder,p5,p11,p12,p6)
    facet!(builder,p6,p12,p7,p1)
    # sides
    facet!(builder,p7,p8,p22,p21)
    facet!(builder,p8,p9,p23,p22)
    facet!(builder,p9,p10,p24,p23)
    facet!(builder,p10,p11,p25,p24)
    facet!(builder,p11,p12,p26,p25)
    facet!(builder,p12,p7,p21,p26)
    # top
    facet!(builder,p15,p21,p22,p16)
    facet!(builder,p16,p22,p23,p17)
    facet!(builder,p17,p23,p24,p18)
    facet!(builder,p18,p24,p25,p19)
    facet!(builder,p19,p25,p26,p20)
    facet!(builder,p20,p26,p21,p15)

    facetregion!(builder,3) # stressor region (bottom, sides, top)
    # bottom
    facet!(builder,p7,p12,p13)
    facet!(builder,p11,p12,p13,p14)
    facet!(builder,p10,p11,p14)
    # sides
    facet!(builder,p7,p13,p27,p21)
    facet!(builder,p13,p14,p28,p27)
    facet!(builder,p10,p14,p28,p24)
    # top
    facet!(builder,p21,p26,p27)
    facet!(builder,p25,p26,p27,p28)
    facet!(builder,p24,p25,p28)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0,0,scale[4]/2))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(scale[1]+scale[2]/2,0,scale[4]/2))

    cellregion!(builder,3) # material 3
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2,scale[4]/2))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 3, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function nanowire_tensorgrid(; scale = [1,1,1,1], nrefs = 1, cut_levels = scale[4]/2, α = nothing, Plotter = nothing, z_levels_dist = 100, version = 1)

    @info "Generating nanowire grid for geometry = $scale"

    builder=SimplexGridBuilder(Generator=Triangulate)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ = scale[3]

    A_core = 3*sqrt(3)/2 * scale[1]^2
    A_shell = 3*sqrt(3)/2 * (scale[2]^2 + 2*scale[1]*scale[2])
    A_stressor = sqrt(3)/2 * scale[3] * (7*(scale[1]+scale[2]) + 3*scale[3])
    if α !== nothing
        A_interface = 3*(d2 * sqrt(3)/2*α)
        A_shell = A_shell - A_interface
        A_stressor = A_stressor - A_interface

        vol_factor_core = 4.0^-1
        vol_factor_shell = 4.0^-1
        vol_factor_interface = 4.0^-(nrefs+1)
        vol_factor_stressor = 4.0^-nrefs
    else
        vol_factor_core = 4.0^-nrefs
        vol_factor_shell = 4.0^-nrefs
        vol_factor_stressor = 4.0^-nrefs
    end
    hz_factor = 2.0^-nrefs

    # bottom side at Z = 0
    # p99= point!(builder,0,-d2)
    p0 = point!(builder,0,0)
    p1 = point!(builder,d1,0)
    p2 = point!(builder,d1/2,sqrt(3)/2*d1)
    p3 = point!(builder,-d1/2,sqrt(3)/2*d1)
    p4 = point!(builder,-d1,0)
    p5 = point!(builder,-d1/2,-sqrt(3)/2*d1)
    p6 = point!(builder,d1/2,-sqrt(3)/2*d1)
    p7 = point!(builder,d2,0)
    p8 = point!(builder,d2/2,sqrt(3)/2*d2)
    p9 = point!(builder,-d2/2,sqrt(3)/2*d2)
    p10 = point!(builder,-d2,0)
    p11 = point!(builder,-d2/2,-sqrt(3)/2*d2)
    p12 = point!(builder,d2/2,-sqrt(3)/2*d2)
    if version == 1
        p13 = point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ)
        p14 = point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ)
    elseif version == 2
        p13 = point!(builder,d2,-δ)
        p14 = point!(builder,d2/2,-sqrt(3)/2*d2-δ)
        p15 = point!(builder,-d2/2,-sqrt(3)/2*d2-δ)
        p16 = point!(builder,-d2,-δ)
    end

    if α !== nothing
        xstar = α/2 + d2*(1 - sqrt(3)*α/(4*δ))
        ystar = - α*(2*sqrt(3)*δ + 3*d2)/(4*δ)
        p15 = point!(builder,-d2+α/2,sqrt(3)/2*α)
        p16 = point!(builder,(α-d2)/2,-sqrt(3)/2*(d2-α))
        p17 = point!(builder,(d2-α)/2,-sqrt(3)/2*(d2-α))
        p18 = point!(builder,d2-α/2,sqrt(3)/2*α)
        p19 = point!(builder,xstar,ystar)
        p20 = point!(builder,(d2+α)/2,-sqrt(3)/2*(d2+α))
        p21 = point!(builder,-(d2+α)/2,-sqrt(3)/2*(d2+α))
        p22 = point!(builder,-xstar,ystar)
    end

    if α !== nothing
        facetregion!(builder,1) # core region
        facet!(builder,p1,p2)
        facet!(builder,p0,p1)
        facet!(builder,p0,p2)

        facet!(builder,p2,p3)
        facet!(builder,p0,p3)

        facet!(builder,p3,p4)
        facet!(builder,p0,p4)

        facet!(builder,p4,p5)
        facet!(builder,p0,p5)

        facet!(builder,p5,p6)
        facet!(builder,p0,p6)

        facet!(builder,p6,p1)

        facetregion!(builder,2) # shell region
        facet!(builder,p15,p16)
        facet!(builder,p16,p17)
        facet!(builder,p17,p18)
        facet!(builder,p18,p8)
        facet!(builder,p8,p9)
        facet!(builder,p9,p15)

        facetregion!(builder,3) # interface region (inside)
        facet!(builder,p15,p10)
        facet!(builder,p10,p11)
        facet!(builder,p11,p12)
        facet!(builder,p12,p7)
        facet!(builder,p7,p18)

        facetregion!(builder,4) # interface region (outside)
        facet!(builder,p7,p19)
        facet!(builder,p19,p20)
        facet!(builder,p20,p21)
        facet!(builder,p21,p22)
        facet!(builder,p22,p10)

        facetregion!(builder,5) # stressor region
        facet!(builder,p19,p13)
        facet!(builder,p13,p14)
        facet!(builder,p14,p22)

        cellregion!(builder,1) # material 1 (core)
        maxvolume!(builder,A_core/12*vol_factor_core)
        regionpoint!(builder,(d1/2,sqrt(3)/4*d1))
        regionpoint!(builder,(0,sqrt(3)/4*d1))
        regionpoint!(builder,(-d1/2,sqrt(3)/4*d1))
        regionpoint!(builder,(d1/2,-sqrt(3)/4*d1))
        regionpoint!(builder,(0,-sqrt(3)/4*d1))
        regionpoint!(builder,(-d1/2,-sqrt(3)/4*d1))

        cellregion!(builder,2) # material 2 (shell)
        maxvolume!(builder,A_shell*vol_factor_shell)
        regionpoint!(builder,(scale[1]+scale[2]/2,0))

        cellregion!(builder,3) # material 3 (inside interface)
        maxvolume!(builder,A_interface*vol_factor_interface)
        regionpoint!(builder,(d2-α/2,0))

        cellregion!(builder,4) # material 4 (outside interface)
        maxvolume!(builder,A_interface*vol_factor_interface)
        regionpoint!(builder,(0,-sqrt(3)/2*(d2+α/2)))

        cellregion!(builder,5) # material 5 (stressor)
        maxvolume!(builder,A_stressor*vol_factor_stressor)
        regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2))
    else
        facetregion!(builder,1) # core region
        facet!(builder,p1,p2)
        facet!(builder,p2,p3)
        facet!(builder,p3,p4)
        facet!(builder,p4,p5)
        facet!(builder,p5,p6)
        facet!(builder,p6,p1)

        facetregion!(builder,2) # shell region
        facet!(builder,p7,p8)
        facet!(builder,p8,p9)
        facet!(builder,p9,p10)
        facet!(builder,p10,p11)
        facet!(builder,p11,p12)
        facet!(builder,p12,p7)

        facetregion!(builder,3) # stressor region
        if version == 1
            facet!(builder,p7,p13)
            facet!(builder,p13,p14)
            facet!(builder,p14,p10)
        elseif version == 2
            facet!(builder,p7,p13)
            facet!(builder,p13,p14)
            facet!(builder,p14,p15)
            facet!(builder,p15,p16)
            facet!(builder,p16,p10)
        end

        cellregion!(builder,1) # material 1
        maxvolume!(builder,A_core/6*vol_factor_core)
        regionpoint!(builder,(0,0))

        cellregion!(builder,2) # material 2
        maxvolume!(builder,A_shell/6*vol_factor_shell)
        regionpoint!(builder,(scale[1]+scale[2]/2,0))

        cellregion!(builder,3) # material 3
        maxvolume!(builder,A_stressor/6*vol_factor_stressor)
        regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2))
    end
    xgrid = simplexgrid(builder)

    hz = z_levels_dist * hz_factor
    if cut_levels !== nothing

        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
        for i = 1 : length(cut_levels)
            index = findfirst(item -> item >= cut_levels[i], z_levels)
            if z_levels_nonuniform[index] == cut_levels[i]
                z_levels_nonuniform[index] = cut_levels[i]-hz/2:hz/4:cut_levels[i]+hz/2
            else
                hz1 = (cut_levels[i]-z_levels_nonuniform[index-1])/2
                hz2 = (z_levels_nonuniform[index]-cut_levels[i])/2
                z_levels_nonuniform[index-1] = z_levels_nonuniform[index-1]:hz1:cut_levels[i]
                z_levels_nonuniform[index]   = cut_levels[i]+hz2:hz2:z_levels_nonuniform[index]
            end
        end
    else
        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
    end
    z_levels_nonuniform = vcat(z_levels_nonuniform...)

    if α !== nothing
        cellregions = xgrid[CellRegions]
        for i = 1 : num_cells(xgrid)
            if cellregions[i] == 3
                cellregions[i] = 2
            end
            if cellregions[i] == 4 || cellregions[i] == 5
                cellregions[i] = 3
            end
        end
        xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset = 5, top_offset = 8)
        # the offsets lead to the following boundary regions:
        # 1 = side core (not seen from outside)
        # 2 = side shell
        # 3 = side interface inside
        # 4 = side interface outside
        # 5 = side stressor
        # 6 = bottom core
        # 7 = bottom shell
        # 8 = bottom stressor
        # 9 = top core
        # 10 = top shell
        # 11 = top stressor
    else
        xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset = 3, top_offset = 6)
        # the offsets lead to the following boundary regions:
        # 1 = side core (not seen from outside)
        # 2 = side shell
        # 3 = side stressor
        # 4 = bottom core
        # 5 = bottom shell
        # 6 = bottom stressor
        # 7 = top core
        # 8 = top shell
        # 9 = top stressor
    end

    return xgrid
end


function bimetal_tensorgrid(; scale = [1,1,1], nrefs = 1, material_border = 0.5)

    @info "Generating bimetal grid for scale = $scale"

    vol_factor = scale[1]*scale[2]*4.0^-nrefs/2
    hz_factor = 2.0^-nrefs

    builder=SimplexGridBuilder(Generator=Triangulate)
    p1=point!(builder,0,0)
    p2=point!(builder,scale[1],0)
    p3=point!(builder,scale[1],scale[2])
    p4=point!(builder,0,scale[2])
    p5=point!(builder,0,material_border*scale[2])
    p6=point!(builder,scale[1],material_border*scale[2])

    facetregion!(builder,10) # interior boundary
    facet!(builder,p5 ,p6)

    facetregion!(builder,1) # outer boundary material 1
    facet!(builder,p5, p1)
    facet!(builder,p1, p2)
    facet!(builder,p2, p6)

    facetregion!(builder,2) # outer boundary material 2
    facet!(builder,p6, p3)
    facet!(builder,p3, p4)
    facet!(builder,p4, p5)

    cellregion!(builder,1)
    maxvolume!(builder,vol_factor)
    regionpoint!(builder,(0.5*scale[1],0.5*material_border*scale[2]))

    cellregion!(builder,2)
    maxvolume!(builder,vol_factor)
    regionpoint!(builder,(0.5*scale[1],0.5*(material_border+1)*scale[2]))

    xgrid = simplexgrid(builder)

    hz = 100 * hz_factor
    xgrid = simplexgrid(xgrid,0:hz:scale[3]; bot_offset = 2, top_offset = 4)

    return xgrid

end


function bimetal_tensorgrid_uniform(; scale = [1,1,1], nrefs = 1, material_border = 0.5, hz = 50)

	W = scale[1]; H = scale[2]; Z = scale[3]
    h1 = W*(1-material_border)
    h2 = W*material_border

    @info "Generating bimetal 3D grid for scale = $scale and middle interface at $material_border (core = $h1, stressor = $h2)"

    # cross-section refinement
    indicator = 1 - abs(2*material_border - 1) # hat function with maximum value equal to 1 and zero at boundaries of [0,1] interval
    nrefs_cross_section = -1*(indicator <= 0.2) + 0*(0.2 < indicator <= 0.6) + 1*(indicator > 0.6)
    if nrefs_cross_section < 0 && nrefs > 0
        nrefs += nrefs_cross_section
    end

    factor = 2.0^-nrefs
    hx = min(h1,h2)*factor
	hy = hx
	hz = min(hz,Z)

    XX = 0:hx:W
    YY = 0:hy:H
    ZZ = Array{Float64,1}(0:hz:Z)

	if nrefs_cross_section >= 0
        # refined points around material_border
        T = Array{Float64,1}(LinRange(h2-hx/2,h2+hx/2,3))

        # Combine the uniform partition and refined points
        XX = sort(unique([XX; T]))
    end

    xgrid = simplexgrid(XX,YY)
    xgrid = uniform_refine(xgrid,nrefs_cross_section)

    # assigning region numbers: core region = 1, stressor region = 2
    cellmask!(xgrid,[h2,0],[h2,H],1)
    cellmask!(xgrid,[0,0],[h2,H],2)

    xgrid_cross_section = deepcopy(xgrid)
    xgrid = simplexgrid(xgrid, ZZ)
    xgrid = uniform_refine(xgrid,nrefs-1)
    # the offsets lead to the following boundary regions:
    # 1 - 6 = sides core & stressor
    # 7 = bottom core
    # 8 = bottom stressor
    # 9 = top core
    # 10 = top stressor

    # boundary faces
    bfacemask!(xgrid,[h2,0,0],[W,0,Z],1)  # side core left
    bfacemask!(xgrid,[W,0,0],[W,H,Z],2)   # side core
    bfacemask!(xgrid,[h2,H,0],[W,H,Z],3)  # side core right
    bfacemask!(xgrid,[0,H,0],[h2,H,Z],4)  # side stressor left
    bfacemask!(xgrid,[0,0,0],[0,H,Z],5)   # side stressor
    bfacemask!(xgrid,[0,0,0],[h2,0,Z],6)  # side stressor right
    bfacemask!(xgrid,[h2,0,0],[W,H,0],7)  # bottom core
    bfacemask!(xgrid,[0,0,0],[h2,H,0],8)  # bottom stressor
    bfacemask!(xgrid,[h2,0,Z],[W,H,Z],9)  # top core
    bfacemask!(xgrid,[0,0,Z],[h2,H,Z],10) # top stressor

    return xgrid, xgrid_cross_section
end


function bimetal_tensorgrid_uniform!(; scale = [1,1,1], nrefs = 1, material_border = 0.5)

    @info "Generating bimetal 3D grid for scale = $scale and middle interface at $material_border of height $(scale[2])"

    W = scale[1]
    H = scale[2]
    Z = scale[3]
    h1 = round(scale[2]*material_border)
    h2 = round(scale[2]*(1 - material_border))
    hz_factor = 2.0^-nrefs

    hx = W/4*2.0^-nrefs
    hy = min(h1,h2)*2.0^-nrefs
    hz = 100 * hz_factor

    XX = 0:hx:W
    #YY = Array{Float64,1}(0:(h1-hy)/2:h1-hy)
    #append!(YY, Array{Float64,1}(LinRange(h1-hy,h1+hy,3)[2:end-1]))
    #append!(YY, Array{Float64,1}(h1+hy:(h2-hy)/2:H))
    YY = 0:hy:H
    ZZ = Array{Float64,1}(0:hz:Z)

    xgrid = simplexgrid(XX,YY)
    xgrid_cross_section = deepcopy(xgrid)
    xgrid = simplexgrid(xgrid, ZZ)
    xgrid = uniform_refine(xgrid,nrefs)
    # the offsets lead to the following boundary regions:
    # 1 - 6 = sides core & bottom
    # 7 = bottom core
    # 8 = bottom stressor
    # 9 = top core
    # 10 = top stressor

    # assigning region numbers: core region = 1, stressor region = 2
    cellmask!(xgrid,[0,0,0],[W,h1,Z],1)
    cellmask!(xgrid,[0,h1,0],[W,H,Z],2)

    # boundary faces
    bfacemask!(xgrid,[0,0,0],[W,h1,Z],1)  # side core left
    bfacemask!(xgrid,[0,h1,0],[0,H,Z],2)  # side stressor left
    bfacemask!(xgrid,[0,H,0],[W,H,Z],3)   # side stressor
    bfacemask!(xgrid,[W,0,0],[W,h1,Z],4)  # side core right
    bfacemask!(xgrid,[W,h1,0],[W,H,Z],5)  # side stressor right
    bfacemask!(xgrid,[0,0,0],[W,0,Z],6)   # side core
    bfacemask!(xgrid,[0,0,0],[W,h1,0],7)  # bottom core
    bfacemask!(xgrid,[0,h1,0],[W,H,0],8)  # bottom stressor
    bfacemask!(xgrid,[0,0,Z],[W,h1,Z],9)  # top core
    bfacemask!(xgrid,[0,h1,Z],[W,H,Z],10) # top stressor

    return xgrid, xgrid_cross_section

end


function nanowire_tensorgrid_mirror(; scale = [1,1,1,1], shape = 1,
    nrefs = 1, z_nrefs = 2, z_levels_dist = 100, cut_levels = scale[4]/2,
    refinement_width = nothing, corner_refinement = false, manual_refinement = false,
    rotate = true, max_nodes = 20)

    @info "Generating nanowire grid for geometry = $scale"

    builder = SimplexGridBuilder(Generator=Triangulate)
    p::Array{Int64,1} = zeros(Int64,max_nodes)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ  = scale[3]

    ## assign nodes
    p = asign_nodes(p,builder,shape,d1,d2,δ)

    ## assign extra nodes for refinement accross the material interface
    if refinement_width !== nothing
        p = interface_refinement(p,builder,shape,d1,d2,δ,refinement_width)

        if manual_refinement == true
            refine(builder,shape,d1,d2,δ,refinement_width)
        end
    end

    ## building edges
    assign_core_shell_edges(p,builder,shape)
    assign_stressor_edges(p,builder,shape,refinement_width)

    ## assigning regions
    assign_regions(p,builder,shape,refinement_width,nrefs,d1,d2,δ)

    # optional refinement at interface corners
    if corner_refinement == true
        function unsuitable(x1,y1,x2,y2,x3,y3, area)
            bary = [(x1+x2+x3)/3,(y2+y2+y3)/3]
            dist = min(norm(bary-refinement_center1),norm(bary-refinement_center2))
            if shape == 3
                dist = min(dist,norm(bary-refinement_center3))
            end
            if area > 1.5*dist
                return 1
            else
                return 0
            end
        end

        if shape == 1
            refinement_center1 = [-sqrt(3)/2*d2,-d2/2]
            refinement_center2 = [0,-d2]
        else
            refinement_center1 = [-d2,0]
            refinement_center2 = [-d2/2,-sqrt(3)/2*d2]
            refinement_center3 = [-d2/2,sqrt(3)/2*d2]
        end
        options!(builder, unsuitable=unsuitable)
    end

    # build grid
    xgrid = simplexgrid(builder)

    # mirror grid
    xgrid_flipped = deepcopy(xgrid)
    xgrid_flipped[Coordinates][1,:] .*= -1
    xgrid = glue(xgrid,xgrid_flipped; interface=4)
    bfacemask!(xgrid,[0,-(d2+δ)],[0,d2+δ],0)
    if shape == 3
        xgrid_flipped = deepcopy(xgrid)
        xgrid_flipped[Coordinates][2,:] .*= -1
        xgrid = glue(xgrid,xgrid_flipped; interface=4)
        bfacemask!(xgrid,[-(d2+2/sqrt(3)*δ),0],[d2+2/sqrt(3)*δ,0],0)
    end

    # clockwise grid rotation by angle θ
	θ = rotate * pi/180
	xgrid_rotated = deepcopy(xgrid)
	xgrid_rotated[Coordinates][1,:] = cos(θ)*xgrid[Coordinates][1,:] +
		sin(θ)*xgrid[Coordinates][2,:]
	xgrid_rotated[Coordinates][2,:] = -sin(θ)*xgrid[Coordinates][1,:] +
		cos(θ)*xgrid[Coordinates][2,:]
	xgrid = deepcopy(xgrid_rotated)
	repair_grid!(xgrid)

    # save cross-section grid
    xgrid_cross_section=deepcopy(xgrid)

    # tensor grid in the z-direction
    hz_factor = 2.0^-z_nrefs
    hz = z_levels_dist * hz_factor
    if cut_levels !== nothing

        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
        for i = 1 : length(cut_levels)
            index = findfirst(item -> item >= cut_levels[i], z_levels)
            if z_levels_nonuniform[index] == cut_levels[i]
                z_levels_nonuniform[index] =
                    cut_levels[i]-hz/2:hz/4:cut_levels[i]+hz/2
            else
                hz1 = (cut_levels[i]-z_levels_nonuniform[index-1])/2
                hz2 = (z_levels_nonuniform[index]-cut_levels[i])/2
                z_levels_nonuniform[index-1] =
                    z_levels_nonuniform[index-1]:hz1:cut_levels[i]
                z_levels_nonuniform[index] =
                    cut_levels[i]+hz2:hz2:z_levels_nonuniform[index]
            end
        end
    else
        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
    end
    z_levels_nonuniform = vcat(z_levels_nonuniform...)

    # generating final grid
    xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=3, top_offset=6)
    # the offsets lead to the following boundary regions:
    # 1 = side core (not seen from outside)
    # 2 = side shell
    # 3 = side stressor
    # 4 = bottom core
    # 5 = bottom shell
    # 6 = bottom stressor
    # 7 = top core
    # 8 = top shell
    # 9 = top stressor

    return xgrid, xgrid_cross_section
end


function asign_nodes(p,builder,shape,d1,d2,δ)
    """
        Assign nodes for different nanowire cross-section shapes.
    """
    if shape == 1
        p[1] = point!(builder,0,0)
        p[2] = point!(builder,0,d1)
        p[3] = point!(builder,-sqrt(3)/2*d1,d1/2)
        p[4] = point!(builder,-sqrt(3)/2*d1,-d1/2)
        p[5] = point!(builder,0,-d1)
        p[6] = point!(builder,0,d2)
        p[7] = point!(builder,-sqrt(3)/2*d2,d2/2)
        p[8] = point!(builder,-sqrt(3)/2*d2,-d2/2)
        p[9] = point!(builder,0,-d2)
        p[10] = point!(builder,0,-d2-δ)
        p[11] = point!(builder,-sqrt(3)/2*d2,-d2/2-δ)
    elseif shape == 2
        p[1] = point!(builder,0,0)
        p[2] = point!(builder,0,sqrt(3)/2*d1)
        p[3] = point!(builder,-d1/2,sqrt(3)/2*d1)
        p[4] = point!(builder,-d1,0)
        p[5] = point!(builder,-d1/2,-sqrt(3)/2*d1)
        p[6] = point!(builder,0,-sqrt(3)/2*d1)
        p[7] = point!(builder,0,sqrt(3)/2*d2)
        p[8] = point!(builder,-d2/2,sqrt(3)/2*d2)
        p[9] = point!(builder,-d2,0)
        p[10] = point!(builder,-d2/2,-sqrt(3)/2*d2)
        p[11] = point!(builder,0,-sqrt(3)/2*d2)
        p[12] = point!(builder,0,-sqrt(3)/2*d2-δ)
        p[13] = point!(builder,-d2/2,-sqrt(3)/2*d2-δ)
        p[14] = point!(builder,-d2,-δ)
    elseif shape == 3
        p[1] = point!(builder,0,0)
        p[2] = point!(builder,0,sqrt(3)/2*d1)
        p[3] = point!(builder,-d1/2,sqrt(3)/2*d1)
        p[4] = point!(builder,-d1,0)
        p[5] = point!(builder,0,sqrt(3)/2*d2)
        p[6] = point!(builder,-d2/2,sqrt(3)/2*d2)
        p[7] = point!(builder,-d2,0)
        p[8] = point!(builder,0,sqrt(3)/2*d2+δ)
        p[9] = point!(builder,-(d2+2/sqrt(3)*δ)/2,sqrt(3)/2*d2+δ)
        p[10] = point!(builder,-d2-2/sqrt(3)*δ,0)
    end

    return p
end


function interface_refinement(p,builder,shape,d1,d2,δ,α)
    """
        Allocate nodes at additonal regions along the material interface.
    """
    if shape == 1
        if α > d2 - d1 || α > δ
            @warn "The refinement width cannot be larger than the shell's diameter or stressor width."
            return nothing
        end

        p[8] = point!(builder,-sqrt(3)/2*d2,-d2/2+α)
        p[9] = point!(builder,0,-d2+α)
        p[10] = point!(builder,0,-d2-δ)
        p[11] = point!(builder,-sqrt(3)/2*d2,-d2/2-δ)
        p[12] = point!(builder,-sqrt(3)/2*d2,-d2/2-α)
        p[13] = point!(builder,0,-d2-α)
        p[14] = point!(builder,0,-d2)
        p[15] = point!(builder,-sqrt(3)/2*d2,-d2/2)
    elseif shape == 2
        if α >= d2 - d1 || α >= δ/sqrt(3)
            @warn "The refinement width cannot be larger than the shell's diameter or δ/sqrt(3)."
            return nothing
        end

        p[9] = point!(builder,-d2+α/2,sqrt(3)/2*α)
        p[10] = point!(builder,(α-d2)/2,-sqrt(3)/2*(d2-α))
        p[11] = point!(builder,0,-sqrt(3)/2*(d2-α))
        p[15] = point!(builder,0,-sqrt(3)/2*d2)
        p[16] = point!(builder,-d2/2,-sqrt(3)/2*d2)
        p[17] = point!(builder,-d2,0)
        p[18] = point!(builder,-d2,-sqrt(3)*α)
        p[19] = point!(builder,-(d2+α)/2,-sqrt(3)/2*(d2+α))
        p[20] = point!(builder,0,-sqrt(3)/2*(d2+α))
    elseif shape == 3
        if α > d2 - d1 || α > δ
            @warn "The refinement width cannot be larger than the shell's or stressor's diameter."
            return nothing
        end
        p[5] = point!(builder,0,sqrt(3)/2*(d2-α))
        p[6] = point!(builder,(α-d2)/2,sqrt(3)/2*(d2-α))
        p[7] = point!(builder,-d2+α,0)
        p[8] = point!(builder,0,sqrt(3)/2*d2)
        p[9] = point!(builder,-d2/2,sqrt(3)/2*d2)
        p[10] = point!(builder,-d2,0)
        p[11] = point!(builder,0,sqrt(3)/2*(d2+α))
        p[12] = point!(builder,-(α+d2)/2,sqrt(3)/2*(d2+α))
        p[13] = point!(builder,-d2-α,0)
        p[14] = point!(builder,0,sqrt(3)/2*d2+δ)
        p[15] = point!(builder,-(d2+2/sqrt(3)*δ)/2,sqrt(3)/2*d2+δ)
        p[16] = point!(builder,-d2-2/sqrt(3)*δ,0)
    end

    return p
end


function refine(builder,shape,d1,d2,δ,α)
    """
        Assign points along the middle of the each interface region to enable a
        uniform refinement.
    """
    if shape == 1
        num_pts = trunc(Int, 4*δ/α)
        for n = 0 : num_pts
            g = n/num_pts
            # convex combination between p8 & p15 midpoint and p9 & p14 midpoint
            px = g*(-sqrt(3)/2*d2)
            py = g*(-d2/2+α/2) + (1-g)*(-d2+α/2)
            point!(builder,px,py)
            # convex combination between p15 & p12 midpoint and p14 & p13 midpoint
            px = g*(-sqrt(3)/2*d2)
            py = g*(-d2/2-α/2) + (1-g)*(-d2-α/2)
            point!(builder,px,py)
        end
    elseif shape == 2
        num_pts = 10 # trunc(Int, 4*δ/(sqrt(3)*α))
        for n = 0 : num_pts
            g = n/num_pts
            ## adding points along the horizontal interface
            # convex combination between p10 & p16 midpoint and p11 & p15 midpoint
            px = g*(α/2-d2)/2
            py = -sqrt(3)/2*(d2-α/2)
            point!(builder,px,py)
            # convex combination between p16 & p19 midpoint and p15 & p20 midpoint
            px = g*(-(d2+α/2)/2)
            py = -sqrt(3)/2*(d2+α/2)
            point!(builder,px,py)
        end
        num_pts = 2*num_pts
        for n = 0 : num_pts-1
            g = n/num_pts
            ## adding points along the left interface
            # convex combination between p9 & p17 midpoint and p10 & p16 midpoint
            px = g*(-d2+α/4) + (1-g)*(α/2-d2)/2
            py = g*sqrt(3)/2*α/2 + (1-g)*(-sqrt(3)/2*(d2-α/2))
            point!(builder,px,py)
            # convex combination between p17 & p18 midpoint and p16 & p19 midpoint
            px = g*(-d2) + (1-g)*(-(d2+α/2)/2)
            py = g*(-sqrt(3)*α/2) + (1-g)*(-sqrt(3)/2*(d2+α/2))
            point!(builder,px,py)
        end
    elseif shape == 3
        num_pts = trunc(Int, d2/4) + 2
        for n = 0 : num_pts
            g = n/num_pts
            # convex combination between p5 & p8 midpoint and p6 & p9 midpoint
            px = (1-g)*(α/2-d2)/2
            py = sqrt(3)/2*(d2-α/2)
            point!(builder,px,py)
            # convex combination between p8 & p11 midpoint and p9 & p12 midpoint
            px = (1-g)*(-α/2-d2)/2
            py = sqrt(3)/2*(d2+α/2)
            point!(builder,px,py)
        end
        num_pts = 2*num_pts
        for n = 0 : num_pts-1
            g = n/num_pts
            # convex combination between p6 & p9 midpoint and p7 & p10 midpoint
            px = g*(α/2-d2)/2 + (1-g)*(-d2+α/2)
            py = g*sqrt(3)/2*(d2-α/2)
            point!(builder,px,py)
            # convex combination between p9 & p12 midpoint and p10 & p13 midpoint
            px = g*(-α/2-d2)/2 + (1-g)*(-d2-α/2)
            py = g*sqrt(3)/2*(d2+α/2)
            point!(builder,px,py)
        end
    end

end


function assign_core_shell_edges(p,builder,shape)
    """
        Assign edges.
    """
    if shape == 1
        facetregion!(builder,1) # core region
        facet!(builder,p[1],p[2])
        facet!(builder,p[2],p[3])
        facet!(builder,p[3],p[4])
        facet!(builder,p[4],p[5])
        facet!(builder,p[5],p[1])

        facetregion!(builder,2) # shell region
        facet!(builder,p[2],p[6])
        facet!(builder,p[6],p[7])
        facet!(builder,p[7],p[8])
        facet!(builder,p[8],p[9])
        facet!(builder,p[9],p[5])
    elseif shape == 2
        facetregion!(builder,1) # core region
        facet!(builder,p[1],p[2])
        facet!(builder,p[2],p[3])
        facet!(builder,p[3],p[4])
        facet!(builder,p[4],p[5])
        facet!(builder,p[5],p[6])
        facet!(builder,p[6],p[1])

        facetregion!(builder,2) # shell region
        facet!(builder,p[2],p[7])
        facet!(builder,p[7],p[8])
        facet!(builder,p[8],p[9])
        facet!(builder,p[9],p[10])
        facet!(builder,p[10],p[11])
        facet!(builder,p[11],p[6])
    elseif shape == 3
        facetregion!(builder,1) # core region
        facet!(builder,p[1],p[2])
        facet!(builder,p[2],p[3])
        facet!(builder,p[3],p[4])
        facet!(builder,p[4],p[1])

        facetregion!(builder,2) # shell region
        facet!(builder,p[2],p[5])
        facet!(builder,p[5],p[6])
        facet!(builder,p[6],p[7])
        facet!(builder,p[7],p[4])
    end
end


function assign_stressor_edges(p,builder,shape,refinement_width)
    """
        Assign stressor edges.
    """
    if refinement_width == nothing
        facetregion!(builder,3) # stressor region
        if shape == 1
            facet!(builder,p[9],p[10])
            facet!(builder,p[10],p[11])
            facet!(builder,p[11],p[8])
        elseif shape == 2
            facet!(builder,p[11],p[12])
            facet!(builder,p[12],p[13])
            facet!(builder,p[13],p[14])
            facet!(builder,p[14],p[9])
        elseif shape == 3
            facet!(builder,p[5],p[8])
            facet!(builder,p[8],p[9])
            facet!(builder,p[9],p[10])
            facet!(builder,p[10],p[7])
        end
    else
        if shape == 1
            facetregion!(builder,2) # interface shell region
            facet!(builder,p[9],p[14])
            facet!(builder,p[14],p[15])
            facet!(builder,p[15],p[8])

            facetregion!(builder,3) # interface stressor region
            facet!(builder,p[15],p[12])
            facet!(builder,p[12],p[13])
            facet!(builder,p[13],p[14])

            facetregion!(builder,3) # stressor region
            facet!(builder,p[13],p[10])
            facet!(builder,p[10],p[11])
            facet!(builder,p[11],p[12])
        elseif shape == 2
            facetregion!(builder,2) # interface shell region
            facet!(builder,p[11],p[15])
            facet!(builder,p[15],p[16])
            facet!(builder,p[16],p[17])
            facet!(builder,p[17],p[9])

            facetregion!(builder,3) # interface stressor region
            facet!(builder,p[17],p[18])
            facet!(builder,p[18],p[19])
            facet!(builder,p[19],p[20])
            facet!(builder,p[20],p[15])

            facetregion!(builder,3) # stressor region
            facet!(builder,p[20],p[12])
            facet!(builder,p[12],p[13])
            facet!(builder,p[13],p[14])
            facet!(builder,p[14],p[18])
        elseif shape == 3
            facetregion!(builder,2) # interface shell region
            facet!(builder,p[5],p[8])
            facet!(builder,p[8],p[9])
            facet!(builder,p[9],p[10])
            facet!(builder,p[10],p[7])

            facetregion!(builder,3) # interface stressor region
            facet!(builder,p[8],p[11])
            facet!(builder,p[11],p[12])
            facet!(builder,p[12],p[13])
            facet!(builder,p[13],p[10])

            facetregion!(builder,3) # stressor region
            facet!(builder,p[11],p[14])
            facet!(builder,p[14],p[15])
            facet!(builder,p[15],p[16])
            facet!(builder,p[16],p[13])
        end
    end
end


function assign_regions(p,builder,shape,refinement_width,nrefs,d1,d2,δ)
    """
        Assign regions.
    """
    A_core = 1/2 * (3*sqrt(3)/2*d1^2)
    A_shell = 1/2 * (3*sqrt(3)/2*(d2^2 - d1^2))

    vol_factor_core = 4.0^-nrefs
    vol_factor_shell = 4.0^-nrefs
    vol_factor_stressor = 4.0^-nrefs

    if shape == 1

        A_stressor = δ*d2

        if refinement_width == nothing
            cellregion!(builder,3) # material 3 (stressor)
            maxvolume!(builder,A_stressor*vol_factor_stressor)
            regionpoint!(builder,(-sqrt(3)/4*d2,-d2))
        else
            α = refinement_width
            A_interface_interior = d2*α
            # the refinement area of the stressor region is the same as the
            # refinement region in the shell
            A_interface_exterior = A_interface_interior
            A_shell = A_shell - A_interface_interior
            A_stressor = A_stressor - A_interface_exterior

            vol_factor_core = 4.0^-(nrefs-1)
            vol_factor_interface = 4.0^-nrefs
            vol_factor_stressor = 4.0^-(nrefs-1)

            cellregion!(builder,2) # material 2 (interface shell)
            maxvolume!(builder,A_interface_interior*vol_factor_interface)
            regionpoint!(builder,(-sqrt(3)/4*d2,1/2*(-3/2*d2+α)))

            cellregion!(builder,3) # material 3 (interface stressor)
            maxvolume!(builder,A_interface_exterior*vol_factor_interface)
            regionpoint!(builder,(-sqrt(3)/4*d2,1/2*(-3/2*d2-α)))

            cellregion!(builder,3) # material 3 (stressor)
            maxvolume!(builder,A_stressor*vol_factor_stressor)
            regionpoint!(builder,(-sqrt(3)/4*d2,-3/4*d2-(δ+α)/2)) # -d2-δ+d2/4+(δ-α)/2
        end

        cellregion!(builder,1) # material 1 (core)
        maxvolume!(builder,A_core*vol_factor_core)
        regionpoint!(builder,(-sqrt(3)/4*d1,0))

        cellregion!(builder,2) # material 2 (shell)
        maxvolume!(builder,A_shell*vol_factor_shell)
        regionpoint!(builder,(-sqrt(3)/2*(d1+d2)/2),0)

    elseif shape == 2

        A_stressor = 1/2 * (δ*(d2 + δ/sqrt(3)) + δ*d2)

        if refinement_width == nothing
            cellregion!(builder,3) # material 3 (stressor)
            maxvolume!(builder,A_stressor*vol_factor_stressor)
            regionpoint!(builder,(-d1/2,-sqrt(3)/2*d2-δ/2))
        else
            α = refinement_width
            A_interface_interior = (2*d2-α)*sqrt(3)/2*α + d2*α
            # the refinement area of the stressor region is almost the same as the
            # refinement region in the shell
            A_interface_exterior = A_interface_interior
            A_shell = A_shell - A_interface_interior
            A_stressor = A_stressor - A_interface_exterior

            vol_factor_core = 4.0^-(nrefs-1)
            vol_factor_interface = 4.0^-nrefs
            vol_factor_stressor = 4.0^-(nrefs-1)

            cellregion!(builder,2) # material 2 (interface shell)
            maxvolume!(builder,A_interface_interior*vol_factor_interface)
            regionpoint!(builder,((α/2-d2)/2,-sqrt(3)/2*(d2-α/2)))

            cellregion!(builder,3) # material 3 (interface stressor)
            maxvolume!(builder,A_interface_exterior*vol_factor_interface)
            regionpoint!(builder,(-(d2+α/2)/2,-sqrt(3)/2*(d2+α/2)))

            cellregion!(builder,3) # material 3 (stressor)
            maxvolume!(builder,A_stressor*vol_factor_stressor)
            regionpoint!(builder,(-d1/2,-sqrt(3)/2*d2-(δ+α)/2))
        end

        cellregion!(builder,1) # material 1 (core)
        maxvolume!(builder,A_core*vol_factor_core)
        regionpoint!(builder,(-d1/2,0))

        cellregion!(builder,2) # material 2 (shell)
        maxvolume!(builder,A_shell*vol_factor_shell)
        regionpoint!(builder,(-d1/2,sqrt(3)/2*(d1+d2)/2))

    elseif shape == 3

        A_core /= 2
        A_shell /= 2
        A_stressor = 1/4 * (3*sqrt(3)/2*((d2+2/sqrt(3)*δ)^2 - d2^2))

        if refinement_width == nothing
            cellregion!(builder,2) # material 2 (shell)
            maxvolume!(builder,A_shell*vol_factor_shell)
            regionpoint!(builder,(-d2/4,sqrt(3)/2*(d1+d2)/2))

            cellregion!(builder,3) # material 3 (stressor)
            maxvolume!(builder,A_stressor*vol_factor_stressor)
            regionpoint!(builder,(-(d2+2/sqrt(3)*δ)/4,sqrt(3)/2*d2+δ/2))
        else
            α = refinement_width
            A_interface_interior = 3/2 * (2*d2-α)*sqrt(3)/2*α
            A_interface_exterior = 3/2 * (2*d2+α)*sqrt(3)/2*α
            A_shell = A_shell - A_interface_interior
            A_stressor = A_stressor - A_interface_exterior

            vol_factor_core = 4.0^-(nrefs-1)
            vol_factor_shell = 4.0^-(nrefs-1)
            vol_factor_interface = 4.0^-nrefs
            vol_factor_stressor = 4.0^-(nrefs-1)

            cellregion!(builder,2) # material 2 (shell)
            maxvolume!(builder,A_shell*vol_factor_shell)
            regionpoint!(builder,(-d2/4,sqrt(3)/2*(d1+d2-α)/2))

            cellregion!(builder,2) # material 2 (interface shell)
            maxvolume!(builder,A_interface_interior*vol_factor_interface)
            regionpoint!(builder,(-d2/4,sqrt(3)/2*(d2-α/2)))

            cellregion!(builder,3) # material 3 (interface stressor)
            maxvolume!(builder,A_interface_exterior*vol_factor_interface)
            regionpoint!(builder,(-d2/4,sqrt(3)/2*(d2+α/2)))

            cellregion!(builder,3) # material 3 (stressor)
            maxvolume!(builder,A_stressor*vol_factor_stressor)
            regionpoint!(builder,(-(d2+2/sqrt(3)*δ)/4,sqrt(3)/2*(2*d2+α)/2+δ/2))
        end

        cellregion!(builder,1) # material 1 (core)
        maxvolume!(builder,A_core*vol_factor_core)
        regionpoint!(builder,(-d1/4,sqrt(3)/2*d1/2))

    end

end
