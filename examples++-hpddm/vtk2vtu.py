#!/usr/bin/env python
# -*- coding: utf-8 -*-

from vtk import *
import os
import sys

base_name = sys.argv[1]
reader = vtkUnstructuredGridReader()
reader.SetFileName(base_name + ".vtk")
reader.Update()
output = reader.GetOutput()

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(base_name + ".vtu")
writer.SetInputData(output)
writer.Write()
os.remove(base_name + ".vtk")
