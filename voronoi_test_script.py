# -*- coding: utf-8 -*-
"""
*****************************************************************************************************
Test the volume method by using the original points rather than the offset points - total volume should sum to 1
*****************************************************************************************************
Created on Mon May 12 09:57:42 2014
The code generates a set of points which represent pores. The voronoi diagram creates
bounded and unbounded cells around the points which are formed by lines that are equidistant from neighbouring points
voronoi_test_4 & 5 assume that the adjoining facet normals point from one pore neighbour to the next within the delaunay
traingulation. This is not the case for a random set of points.

The method for identifying a facet is as follows
For each point
    Check if contained within bounded region
    if yes
        Find neighbours
        for each neighbour get vertices for surrounding region (whether bounded or not)
        identify shared vertices between original and neighbour point - these will form the facet
        These will be orientated in 3d space randomly (and may or may not form a unique facet - check this)
        # Replaced with a rotation # identify the spread of values in each direction and select the axis with smallest spread to delete
        compute convex hull in 2d which will give the connections between vertices to form the facet
        draw connections in 3d space
    in no
        ignore
    include some logic to not repeat neighbour connections i.e store the 2 neighbours in a neighbour array
Repeat until pores in bounded regions have had all facets shared with neighbours identified
@author: pmtgt
"""
" ********************************************************************************************************************************** "    

def centroid(points):
    import numpy as np
    columns = np.shape(points)[1]
    if columns == 2:
        cen = np.array([points[:,0].mean(),points[:,1].mean()])
    elif columns == 3:
        cen = np.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
    else:
        cen = 0
    return cen
        

def prob_func(m):
    a = 35
    b = 0.1
    p = ((m**a) + ((1-m)**a) + (2*b))/(1 + (2*b))
    #p = ((m**a) + b)/(1 + b)
    return p

def reject(point):
    
    import numpy as np
    x = point[0]
    y = point[1]
    z = point[2]
    Px = prob_func(x)
    Py = prob_func(y)
    Pz = prob_func(z)
    nrand = np.random.uniform(0,1,1)
    
    if Px < nrand and Py < nrand:
        rejection = True
    else:
        rejection = False
        #print x, Px, nrand
    
    return rejection

def pore_volume(points,plot=False):
    import numpy as np
    from scipy.spatial import Delaunay

    #points = np.random.randn(20, 3) # Random points in space
    #points = np.array([[0,0,0],[0,1,0],[1,0,0],[1,1,0],[0.5,0.5,3]]) # Regular Conical Pyramid - Volume = 1.0
    tri = Delaunay(points)
    " We only want points included in the convex hull to calculate the centroid "
    hull_points = np.unique(tri.convex_hull)
    #hull_centroid = np.array([points[hull_points,0].mean(),points[hull_points,1].mean(),points[hull_points,2].mean()])
    hull_centroid = centroid(points[hull_points])
    # -- Make a list of faces, [(p1, p2, p3), ...];  pj = (xj, yj, zj)

    faces = [] # list of points making each face
    face_centers = [] # co-ordinates of the face centroids
    face_normals = [] # normal vector of the faces
    face_areas = [] # Area of each face
    hull_volume = 0.0
    for ia, ib, ic in tri.convex_hull:
        " Points making each triangular face "
        faces.append(points[[ia, ib, ic]])
        " Collection of co-ordinates of each point in this face "
        face_x = points[[ia,ib,ic]][:,0]
        face_y = points[[ia,ib,ic]][:,1]
        face_z = points[[ia,ib,ic]][:,2]
        " Average of each co-ordinate is the centroid of the face "
        face_centroid = [face_x.mean(),face_y.mean(),face_z.mean()]
        face_centers.append(face_centroid)
        face_centroid_vector = face_centroid - hull_centroid
        " Vectors of the sides of the face used to find normal vector and area "
        vab = points[ib] - points[ia]
        vac = points[ic] - points[ia]
        vbc = points[ic] - points[ib] # used later for area
        " As vectors are co-planar the cross-product will produce the normal vector of the face "
        face_normal = np.cross(vab,vac)
        face_unit_normal = face_normal/np.linalg.norm(face_normal)
        face_normals.append(face_unit_normal)
        " As triangles are orientated randomly in 3D we could either transform co-ordinates to align with a plane and perform 2D operations "
        " to work out the area or we could work out the lengths of each side and use Heron's formula which is easier"
        " Using Delaunay traingulation will always produce triangular faces but if dealing with other polygons co-ordinate transfer may be necessary "
        a = np.linalg.norm(vab)
        b = np.linalg.norm(vbc)
        c = np.linalg.norm(vac)
        " Semiperimeter "
        s = 0.5*(a+b+c)
        face_area = np.sqrt(s*(s-a)*(s-b)*(s-c))
        face_areas.append(face_area)
        " Now the volume of the pyramid section defined by the 3 face points and the hull centroid can be calculated "
        pyramid_volume = np.abs(np.dot(face_centroid_vector,face_unit_normal)*face_area/3)
        " Each pyramid is summed together to calculate the total volume "
        hull_volume += pyramid_volume
    
    face_centers = np.asarray(face_centers)
    face_normals = np.asarray(face_normals)
    face_areas = np.asarray(face_areas)

    if (plot):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        r,g,b=np.random.rand(3,1)
        items = Poly3DCollection(faces, facecolors=[(r, g, b, 0.1)])
        ax.add_collection(items)
        #ax.scatter(points[hull_points,0], points[hull_points,1], points[hull_points,2], 'o')
        ax.scatter(points[:,0], points[:,1], points[:,2], 'o')
        ax.scatter(hull_centroid[0],hull_centroid[1],hull_centroid[2],c='r',marker='o')
        #ax.scatter(face_centers[:,0], face_centers[:,1], face_centers[:,2],c='g',marker='o')
        plt.show()
        
    print "Pore Volume: "+str(hull_volume)
    return hull_volume
" ********************************************************************************************************************************** "
def PolyArea2D(pts):
    lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
    area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
    return area
" ********************************************************************************************************************************** "
def dist2(p1, p2):
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2
" ********************************************************************************************************************************** "
def fuse(points, d):
    ret = []
    d2 = d * d
    n = len(points)
    taken = [False] * n
    for i in range(n):
        if not taken[i]:
            count = 1
            point = [points[i][0], points[i][1]]
            taken[i] = True
            for j in range(i+1, n):
                if dist2(points[i], points[j]) < d2:
                    point[0] += points[j][0]
                    point[1] += points[j][1]
                    count+=1
                    taken[j] = True
            point[0] /= count
            point[1] /= count
            ret.append((point[0], point[1]))
    return ret
" ********************************************************************************************************************************** "
def offset_vertex(points,rad = 0.01):
    import numpy as np
    import _transformations as tr
    from math import atan2
    debug=False 
    " We are passed in a set of 3 points forming vertices of two adjoining simplexes of the convex hull of a voronoi facet "
    " We need to offset the vertices normal to the fibre direction (or adjoining vectors) by the fibre radius "
    " This is achieved by finding the half angle between the two adjoining vectors and a direction "
    " Mid-point must be the first in the array "
    p0 = np.array(points[0])
    p1 = np.array(points[1])
    p2 = np.array(points[2])
    " Now make the midpoint the origin "
    vector1 = p1-p0
    vector2 = p2-p0

    " Find what quadrant the vector is pointing in - atan2 function takes account of signs "
    " 0 means aligned with x-axis, pi is aligned with -xaxis, positive numbers are positive y and negative numbers are negative y "
    " The angle between the vectors should always be within 180 degrees of one another in a convex hull "

    q1 = atan2(vector1[1],vector1[0])
    q2 = atan2(vector2[1],vector2[0])
    alpha = 0.5*tr.angle_between_vectors(vector1,vector2)
    if (debug):
        if (q1 >=0) & (q1 <= 0.5*np.pi):
            "We are in 1st quadrant"
            print "Vector 1 is in first quandrant"
        elif (q1 >=0) & (q1 <= np.pi):
            "We are in 2nd quadrant"
            print "Vector 1 is in second quandrant"
        elif (q1 >= -np.pi) & (q1 <= -0.5*np.pi):
            "We are in third quadrant"
            print "Vector 1 is in third quandrant"
        else:
            "We are in 4th quadrant"
            print "Vector 1 is in fourth quandrant"
    
        if (q2 >=0) & (q2 <=0.5*np.pi):
            "We are in 1st quadrant"
            print "Vector 2 is in first quandrant"
        elif (q2 >=0) & (q2 <=np.pi):
            "We are in 2nd quadrant"
            print "Vector 2 is in second quandrant"
        elif (q2 >= -np.pi) & (q2 <= -0.5*np.pi):
            "We are in third quadrant"
            print "Vector 2 is in third quandrant"
        else:
            "We are in 4th quadrant"
            print "Vector 2 is in fourth quandrant"
        print "Vector 1: " +str(vector1)
        print "Vector 2: " +str(vector2)
        print "Half angle alpha: "+str(np.around(180*alpha/np.pi,1)) +" degrees"    
        print "Theta 1: "+str(np.around(180*q1/np.pi,1)) +" degrees" 
        print "Theta 2: "+str(np.around(180*q2/np.pi,1)) +" degrees"
    
    " We always want to offset from the first vertex we get to - going anti-clockwise from the x-axis "
    " Check if both vectors point up or both point down - if so the first one we get to will have smaller q value "
    if q1*q2 >=0.0:
        if q1<q2:
            theta = q1
        else:
            theta = q2
    else:
        "if vector 1 is more rotated away from positive xaxis than vector 2 is rotated away from negative xaxis - use it"
        " and vice-versa "
        if (abs(q1)+abs(q2)>np.pi):
            "vectors are pointing negative x so take whichever has positive q-value - like a pacman facing left"
            if q1>=0:
                theta = q1
            else:
                theta = q2
        else:
            "vectors are pointing positive x so take whichever is negative"
            if q1<=0:
                theta = q1
            else:
                theta = q2
    
    x = rad*np.cos(alpha+theta)/np.sin(alpha)
    y = rad*np.sin(alpha+theta)/np.sin(alpha)

    "Add the midpoint back in"
    output = [x+p0[0],y+p0[1]]    
    
    return output
" ********************************************************************************************************************************** "
def rotate_and_print(facet, normal, axis, precision=10,plot=False,fibre_rad=0.02):
    import _transformations as tr
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.spatial import ConvexHull

    throat_area = 0.0
    throat_centroid = []
    offset_verts_3D = []
    output_offset = []
    " For boundaries some facets will already be aligned with the axis - if this is the case a rotation is unnecessary and could also cause problems "
    angle = tr.angle_between_vectors(normal,axis)
    if (angle==0.0)or(angle==np.pi):
        "We are already aligned"
        output = np.around(facet,precision)
        rotate_input = False
    else:
        rotate_input = True
        M = tr.rotation_matrix(tr.angle_between_vectors(normal,axis),tr.vector_product(normal,axis))
        rotated_facet = np.dot(facet,M[:3,:3].T)
        #output = np.around(tr.unit_vector(rotated_facet),precision)
        output = np.around(rotated_facet,precision)
    

    x = output[:,0]
    y = output[:,1]
    z = output[:,2]
    if (np.around(z.std(),3)!=0.000):
        print "Rotation failed"
    facet_coords_2D = np.column_stack((x,y)) ##THIS NEEDS REVISING IF WE WANT TO GENERALISE TO ANY AXIS BUT IT DOESN'T REALLY MATTER
    hull = ConvexHull(facet_coords_2D)
    
    " Work out span of points and set axes scales to cover this and be equal in both dimensions "
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    if (x_range > y_range):
        my_range = x_range
    else:
        my_range = y_range
           
    lower_bound_x = x.min() - my_range*0.1
    upper_bound_x = x.min() + my_range*1.1
    lower_bound_y = y.min() - my_range*0.1
    upper_bound_y = y.min() + my_range*1.1

    " Now we want to effectively erode the facet the fibre radius to simulate the fibre" 
    " We need to check whether any vertices lie very close together and merge them if they do otherwise offsetting will not work "
    " Also if the range in values of the facet is less than the fibre diameter the facet is two small and should be ignored "
    tolerance = my_range*0.1
    verts_2D = facet_coords_2D[hull.vertices]
    fused_verts = fuse(verts_2D,tolerance)
    if (len(fused_verts) <3):
        print "Error: Fused Too Many Verts"
    elif(my_range < fibre_rad*2):
        print "Error: Facet Too small to Erode"
    else:
        offset = []
        for i,vert in enumerate(fused_verts):
            " Collect three adjacent points and compute the offset of the first "
            triplet = (vert, np.roll(fused_verts,-1,axis=0)[i],np.roll(fused_verts,1,axis=0)[i])
            offset.append(offset_vertex(triplet,fibre_rad))
        offset = np.asarray(offset)
        " At this point if everything went well the offset points should lie within the original convex hull "
        " If they don't then the offset value may be greater than the spread of points or we may still "
        " have some points very close together "
        " Make another convex hull including offset points and see whether the area has increased "
        " If this happens the throat is fully eroded and will have zero area "
            
        original_area = PolyArea2D(verts_2D)            
        all_points = np.concatenate((verts_2D,offset),axis=0)
        total_hull = ConvexHull(all_points)            
        total_area = PolyArea2D(all_points[total_hull.vertices])
        if (total_area>original_area):
            print "Error: Offset area larger than original"
        else:
            offset_hull = ConvexHull(offset)
            offset_verts_2D = offset[offset_hull.vertices]
            throat_area = PolyArea2D(offset_verts_2D)
            print "Throat Area: "+str(throat_area)
            " Make 3D again in rotated plane "
            offset_verts_3D = np.column_stack((offset_verts_2D,z[0:len(offset_verts_2D)]))
            " Get matrix to un-rotate the co-ordinates back to the original orientation if we rotated in the first place"
            if (rotate_input):
                M1 = tr.inverse_matrix(M)
                " Unrotate the offset coordinates "
                output_offset = np.dot(offset_verts_3D,M1[:3,:3].T)
            else:
                output_offset = offset_verts_3D
            " Calculate the centroid of the facet to calculate pore to pore distance later "
            throat_centroid = centroid(output_offset)
    
    if (plot):
        temp_fig = plt.figure()
        plt.axis((lower_bound_x,upper_bound_x,lower_bound_y,upper_bound_y))
        " Plot the convex Hull of the original points "
        for simplex in hull.simplices:
            plt.plot(output[simplex,0], output[simplex,1], 'k-',linewidth=2)
        plt.scatter(facet_coords_2D[:,0], facet_coords_2D[:,1])
        " Plot the convex Hull of the offset points "
        if (throat_area>0.0):       
            for offset_simplex in offset_hull.simplices:
                plt.plot(offset[offset_simplex,0], offset[offset_simplex,1], 'g-',linewidth=2)
            plt.scatter(offset[:,0], offset[:,1])

        temp_fig.show()
        
    return output,output_offset,throat_area,throat_centroid
" ********************************************************************************************************************************** "
" Start of main code "
" ********************************************************************************************************************************** "
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull
#from scipy.spatial import Delaunay
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.special import cbrt
import time
import csv
#import _transformations as tr
start_time = time.time()

" This sets the number of pores in the domain - random numbers are between 0 and 1 and the network is approximately cubic but with higher density in the centre "
" Different distributions could be used here and scaling the domain in one direction could provide anisotropy "
" As the length scales are normalised to a 1,1,1 cube we need to work out the real size of the domain to calculate the correct fibre radius "
" Real pores are of diameter 10 - 15 microns and fibre diameter are 7 - 8 microns so about half the pore diameter "
n=100
" Adjust point distribution to seed edges with slightly more pores - to reduce edge volumes when reflected "
mylist = []
iteration = 0
while len(mylist)<n:
    iteration += 1
    point = np.random.rand(3,1)  
    if reject(point) == False:
        mylist.append(point)
mylist = np.array(mylist)
points = np.hstack((np.vstack(mylist[:,0]),np.vstack(mylist[:,1]) ,np.vstack(mylist[:,2])))      
print iteration


domain_length = 1
#points = np.random.rand(n,3)*domain_length
pore_diameter = 2*domain_length/(3*cbrt(n)+1)
fibre_radius = pore_diameter/3 #factor of 10 necessary to get decent porosity

xaxis = np.array([1, 0, 0])
yaxis = np.array([0, 1, 0])
zaxis = np.array([0, 0, 1])

do_reflection = True
if do_reflection:
    " Reflect in X = 1 and 0 "
    Pxp = points.copy()
    Pxp[:,0]=(2-Pxp[:,0])
    Pxm= points.copy()
    Pxm[:,0] = Pxm[:,0]*(-1)
    " Reflect in Y = 1 and 0 "
    Pyp = points.copy()
    Pyp[:,1]=(2-Pxp[:,1])
    Pym = points.copy()
    Pym[:,1] = Pxm[:,1]*(-1)
    " Reflect in Z = 1 and 0 "
    Pzp = points.copy()
    Pzp[:,2]=(2-Pxp[:,2])
    Pzm = points.copy()
    Pzm[:,2] = Pxm[:,2]*(-1)
    all_points = np.vstack((points,Pxp,Pxm,Pyp,Pym,Pzp,Pzm))
else:
    all_points = points.copy()
    
vor = Voronoi(all_points)

faces = []
indexed_points = enumerate(points)
facet_count = 0

unbounded = set([-1])
bounded_sets = []
done_facets = []
" Get the region indexes for our original set of points "
original_regions = vor.point_region[0:n]

total_pore_volume = 0.0
total_domain_volume = 0.0
"Pore Stats Lists"
coordination_number = []
collect_pore_vols = []
collect_throat_areas = []
collect_pore2pore = []
"Setup csv for output "
dataCSV = open('Voronoi_Network_Data.csv', 'w') 
writer = csv.writer(dataCSV, dialect='excel')
writer.writerow(['Data','Pore', 'Throat', 'Neighbour' , 'Throat Area', 'Pore Volume', 'X', 'Y','Z']) 

for i,my_region in enumerate(original_regions):
    polygon = vor.regions[my_region]    
    polygon_set = set(polygon)
    pore_vol = 0.0
    polygon_volume = 0.0
    " if points were reflected to create bound square domain then all should be bounded, if not then we only "
    " want to deal with bounded regions "
    if (len(polygon_set.intersection(unbounded))==0)&(len(polygon_set)>0):
        " Loop through all regions to find neighbours to my_region identified by sharing vertices "
        throat_count = 0
        occluded_throat = 0
        pore_offset_coords = []        
        for j,neighbour_region in enumerate(vor.point_region):            
            if neighbour_region != my_region:
                neighbour_poly = vor.regions[neighbour_region] 
                if all(neighbour_poly != polygon):         
                    neighbour_set = set(neighbour_poly)
                    shared_verts = np.array(list(polygon_set.intersection(neighbour_set)))
                    " We need at least 3 shared vertices to make a facet "
                    if (len(shared_verts)>=3):
                        done_facet = False
                        " Don't want to be printing facets twice "
                        for facet in done_facets:
                            if all(facet==shared_verts):
                                done_facet = True 

                        throat_count +=1
                        print "Pore: "+str(i)+" Throat: "+str(throat_count) + " Connected to Pore: "+str(j)
                        " To print or work with throat facet we need to eliminate a co-oordinate by rotating "
                        coords = np.around(vor.vertices[shared_verts],10)
                        #normal_vector = points[j] - points[i]
                        normal_vector = all_points[j] - all_points[i]
                        pore2pore = np.linalg.norm(normal_vector)                 
                        unit_normal = normal_vector/pore2pore
                        rotated_coords,offset_coords,throat_area,throat_centroid = rotate_and_print(coords,normal_vector,zaxis,10,plot=False,fibre_rad = fibre_radius) 
                        " The rotated_coords should all be coplanar in the zaxis so we can remove this co-ordinate to make 2D "
                        c2c = 0.0
                        if throat_area == 0.0:
                            occluded_throat += 1
                        elif (done_facet == False):
                            collect_throat_areas.append(throat_area)
                            " Calculate the length of facet from pore centroid to neighbour pore centroid "
                            c2c = np.linalg.norm(all_points[j]-throat_centroid)+np.linalg.norm(all_points[i]-throat_centroid)
                            #collect_pore2pore.append(pore2pore)
                            collect_pore2pore.append(c2c)
                        writer.writerow(['Throat',i, throat_count, j, throat_area,c2c])    
                        x = rotated_coords[:,0]
                        y = rotated_coords[:,1]
                        z = rotated_coords[:,2]
                        if (np.around(z.std(),3)!=0.000):
                            print "Rotation failed"
                        facet_coords_2D = np.column_stack((x,y))
                        hull = ConvexHull(facet_coords_2D)
                        if not done_facet:
                            faces.append(coords[hull.vertices])
                            facet_count += 1
                            done_facets.append(shared_verts)
                        if len(offset_coords) >=3:
                            pore_offset_coords.append(offset_coords)
                            
                        #done_facets.append(shared_verts)
        
        " Plot Delaunay Triangulation of all offset coords to replicate pore volume "        
        if len(pore_offset_coords) >=3.0: # we need at least 3 facets to make a volume
            pore_points=[]
            "We now have a list of lists and need one array with all the points "
            for face in pore_offset_coords:
                for vert in face:
                    pore_points.append(vert)
            pore_points = np.asarray(pore_points)        
            pore_vol = pore_volume(pore_points,plot=False)           
            total_pore_volume += pore_vol
        " Now check the volume of the region "
        if len(polygon)>=4: #we need at least 4 vertices to make a volume
            polygon_volume = pore_volume(vor.vertices[polygon], plot=False)
            total_domain_volume += polygon_volume
        if (pore_vol>polygon_volume):
            print "Verts of reduced hull:"            
            print pore_points
            print "Reduced Volume: "+str(pore_vol)
            print "Original Verts:"
            print vor.vertices[polygon]
            print "Original Volume: "+str(polygon_volume)
            #sys.exit("Volume of reduced hull is larger than original")
        print "**************************************************"
        writer.writerow(['Pore',i, throat_count-occluded_throat, 0.0, 0.0,pore_vol, points[i][0],points[i][1],points[i][2]])
        coordination_number.append(throat_count-occluded_throat)
        collect_pore_vols.append(pore_vol)
        
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(points[:,0], points[:,1], points[:,2], 'r.')
items = Poly3DCollection(faces,linewidths=1, alpha=0.2)
face_colours=[(0.5, 0, 0.5, 0.05)]
items.set_facecolor(face_colours)
ax.add_collection(items)
#ax.scatter(vor.vertices[:,0], vor.vertices[:,1], vor.vertices[:,2], 'o')
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_zlim(0,1)
plt.show()

dataCSV.close()

" Plot histograms of pore stats"
plt.figure()
plt.title("Histogram of Pore Volumes")
num,bins,patches = plt.hist(collect_pore_vols,bins=20,normed=1,histtype='bar',color='r')
plt.show()

plt.figure()
plt.title("Histogram of Throat Areas")
num,bins,patches = plt.hist(collect_throat_areas,bins=20,normed=1,histtype='bar',color='b')
plt.show()

plt.figure()
plt.title("Histogram of Centre to Centre Distances")
num,bins,patches = plt.hist(collect_pore2pore,bins=20,normed=1,histtype='bar',color='y')
plt.show()

plt.figure()
plt.title("Histogram of Pore Coordination Number")
num,bins,patches = plt.hist(coordination_number,bins=len(np.unique(coordination_number)),normed=1,histtype='bar',color='g')
plt.show()

print "Number of Points: " +str(n)
print "Number of Throats: " +str(facet_count)
print "Fibre Radius: "+str(fibre_radius)
print "Total Pore Volume: "+str(total_pore_volume)
print "Total Domain Volume: "+str(total_domain_volume)
print "Porosity: " +str(np.around((total_pore_volume/total_domain_volume),3)*100)+"%"
print "Simulation time: " +str(time.time()-start_time)


