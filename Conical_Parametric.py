# -*- coding: utf-8 -*-
"""
@author: Christoph Maier, c.maier@tugraz.at
last changes: 22.04.2024

Script for different parametric conical inductor scripts.
They are not perfect, but are a starting point for further investigations.

The .stl export and the methods to make it watertight can be found here:
https://stackoverflow.com/questions/67963814/how-to-save-a-vtk-render-as-a-vtk-or-stl

"""


import numpy as np
import matplotlib.pyplot as plt
import vtk
import csv
from stl import mesh


"""

The function takes the starting radius, stopping radius (0 for a pointy inductor),
the number of turns and the shape of the inducotr ('exponential' or 'conical').

The height is determined by the parameters.

"""
def inductor_equidistant(starting_radius, stopping_radius, num_turns, points_per_turn, exp_factor, shape):
    total_points = num_turns * points_per_turn
    t = np.linspace(0, num_turns * 2 * np.pi, total_points)

    if shape == 'conical':
        # Linear interpolation between starting and stopping radii
        radius = np.linspace(starting_radius, stopping_radius, total_points)
    elif shape == 'exponential':
        # Exponential interpolation between starting and stopping radii
        b = np.log(stopping_radius / (starting_radius/exp_factor)) / (num_turns * 2 * np.pi)
        radius = starting_radius * np.exp(b * t)
    else:
        raise ValueError("Shape must be 'conical' or 'exponential'")

    # Calculate x and y coordinates
    x = radius * np.cos(t)
    y = radius * np.sin(t)
    z = np.linspace(0, num_turns, total_points)  # Uniform z distribution along the turns
    
    return x,y,z
    
    
    
"""
The function takes the starting radius, stopping radius (0 for a pointy inductor),
the number of turns, its height and two factors for the shape.

The distance between the windings is determined by the two factors. Also the 
shape of the inductor.

"""
def inductor_fixed_height(starting_radius, stopping_radius, num_turns, height, points_per_turn, exp_factor, decreasing_factor, shape='conical'):
    total_points = num_turns * points_per_turn
    t = np.linspace(0, num_turns * 2 * np.pi, total_points)

    # Exponential interpolation between starting and stopping radii
    b = np.log(stopping_radius / (starting_radius/exp_factor)) / (num_turns * 2 * np.pi)
    radius = starting_radius * np.exp(b * t)

    # Calculate x and y coordinates
    x = radius * np.cos(t)
    y = radius * np.sin(t)
    z = spaced_linspace(0, height, total_points, decreasing_factor)

    return x,y,z
    

"""
Generate a vector of 'num' evenly spaced numbers between 'start' and 'stop',
with decreasing or increasing spacing between the numbers, ensuring the
sequence ends at 'stop'.

Parameters:
    start (float): The starting value of the sequence.
    stop (float): The end value of the sequence.
    num (int): Number of samples to generate.
    decrease_factor (float): Factor by which the spacing decreases between numbers.

Returns:
    numpy.ndarray: An array containing 'num' evenly spaced numbers.
"""
def spaced_linspace(start, stop, num, decrease_factor):
    if num == 1:
        return np.array([start])

    # Calculate initial distances based on a reverse cumulative increase
    total_length = stop - start
    weights = np.geomspace(1, decrease_factor, num - 1)[::-1]
    normalized_weights = weights / weights.sum()
    distances = normalized_weights * total_length

    # Calculate the new positions by summing the scaled distances
    new_positions = np.cumsum(np.insert(distances, 0, start))
    return new_positions


"""
TO-DO: Add description of the .stl function
"""
def create_tube_along_curve(points, radius, num_sides, filename):
    # Create an empty list to store the mesh triangles
    tube_faces = []

    # Helper function to create a circle of points in the plane perpendicular to the tangent vector at the curve
    def circle_points(center, normal, radius, num_points):
        # Find two vectors orthogonal to the normal
        if np.allclose(normal, [0, 0, 1]) or np.allclose(normal, [0, 0, -1]):
            v = np.array([1, 0, 0])
        else:
            v = np.array([0, 0, 1])
        u = np.cross(normal, v)
        v = np.cross(u, normal)
        u /= np.linalg.norm(u)
        v /= np.linalg.norm(v)
        
        # Generate points in the circle
        t = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        circle_pts = center[:, np.newaxis] + radius * (np.cos(t) * u[:, np.newaxis] + np.sin(t) * v[:, np.newaxis])
        return circle_pts.T

    # Generate the tube
    for i in range(1, len(points)):
        p0 = points[i - 1]
        p1 = points[i]
        tangent = p1 - p0
        tangent /= np.linalg.norm(tangent)  # Normalize the tangent vector
        
        circle0 = circle_points(p0, tangent, radius, num_sides)
        circle1 = circle_points(p1, tangent, radius, num_sides)
        
        # Create quads between the points
        for j in range(num_sides):
            next_index = (j + 1) % num_sides
            # Each quad is made of two triangles
            triangle1 = [circle0[j], circle0[next_index], circle1[j]]
            triangle2 = [circle1[j], circle0[next_index], circle1[next_index]]
            tube_faces.append(triangle1)
            tube_faces.append(triangle2)

    # Create the mesh object
    tube_mesh = mesh.Mesh(np.zeros(len(tube_faces), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(tube_faces):
        tube_mesh.vectors[i] = np.array(f)
        
    tube_mesh.save(filename)
    print(f"Start Smoothing of STL file '{filename}'.")
    
    #make mesh watertight
    # read in .stl file to 
    mesh_input = vtk.vtkSTLReader()
    mesh_input.SetFileName(filename)
    mesh_input.Update()
    mesh_input = mesh_input.GetOutput()
    
    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputData(mesh_input)
    appendFilter.Update()       
        
    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
    cleanFilter.Update()
    
    # newData = cleanFilter
    fill = vtk.vtkFillHolesFilter()
    fill.SetInputConnection(appendFilter.GetOutputPort())   
    fill.SetHoleSize(100)
    fill.Update()
    
    # new .stl output
    stl_output = vtk.vtkSTLWriter()
    stl_output.SetInputConnection(fill.GetOutputPort())
    stl_output.SetFileName(filename)
    stl_output.Update()
    
    print(f"STL file '{filename}' has been created.")

    



"""
Main function. For test purposes.
"""
if __name__ =='__main__':

    """
    Example useage of the equidistant function.
    
    Parameters:
        starting_radius: Starting radius of the inductor.
        stopping_radius: Stopping radius of the inductor (0.5 for a more or less pointy inductor).
        num_turns: Number of turns of the inductor.
        smoothness: Smoothing factor (defines the number of points for the vectors).
        shape: 'exponential' or 'conical'. Defines the shape of the inductor.
        exponential_factor: For the exponential inductor only. Can control the shape of the inductor
    """
    starting_radius = 50
    stopping_radius = 5
    num_turns = 10
    smoothness =1000
    exponential_factor = 1 # 1 is standard
    shape = 'exponential'
    #shape = 'conical'
    
    [x1, y1, z1] = inductor_equidistant(starting_radius, stopping_radius, num_turns, smoothness, exponential_factor, shape)
    
    
    """
    Example useage of the fixed height exponential function.
    
    Parameters:
        starting_radius: Starting radius of the inductor.
        stopping_radius: Stopping radius of the inductor (0.5 for a more or less pointy inductor).
        num_turns: Number of turns of the inductor.
        height: Height of the inductor.
        smoothness: Smoothing factor (defines the number of points for the vectors).
        exponential_factor: For the exponential inductor only. Can control the shape of the inductor
        decreasing_factor: Defines how much the gradient of the inductors should change (1 for equidistant, >1 for increasing winding to the tip, >1 for a decreasing winding )
        
        exponential_factor and decreasing_factor define the shape of the inductor. To get a desired shape a little bit of trying is needed. I was not
        able to define that one on a better way.
    """
    starting_radius = 50
    stopping_radius = 5
    num_turns = 10
    height = 70
    smoothness = 1000
    exponential_factor = 0.6 # 1 is standard
    decreasing_factor = 3
    
    [x2, y2, z2] = inductor_fixed_height(starting_radius, stopping_radius, num_turns, height, smoothness, exponential_factor, decreasing_factor, shape)
    
    
    """
        3D and 2D Plotting of the results
    """
    # 3D-plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x1, y1, z1)
    ax.plot(x2, y2, z2)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.title('Different conical forms')
    plt.show()
    
    # 2D-plot
    fig = plt.figure()
    plt.plot(x1, z1, label='equidistant')
    plt.plot(x2, z2, label='variable distance')
    plt.grid()
    plt.legend()
    
    """
        generate vector text file and export it
    """
    
    with open('Coordinates_equidistant.csv', 'w', newline='') as csvfile:
        fieldnames = ['x-value', 'y-value','z-value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=';')
    
        writer.writeheader()
        for i in range(len(x1)):
        
            writer.writerow({'x-value': x1[i], 'y-value': y1[i], 'z-value': z1[i]})
    
    
    
    
    """
        generate .stl file and export it
    """
    points = np.vstack((x1, y1, z1)).T
    radius = 1
    smoothness_stl = 2
    filename='curve_tube_equidistant.stl'
    
    create_tube_along_curve(points, radius, smoothness_stl, filename)
    
    
    
    points = np.vstack((x2, y2, z2)).T
    radius = 1
    smoothness_stl = 2
    filename='curve_tube_variable.stl'
    
    create_tube_along_curve(points, radius, smoothness_stl, filename)
    
    
    
    
    
    