import sys
import vtk

# https://visit-sphinx-github-user-manual.readthedocs.io/en/3.4rc/data_into_visit/VTKFormat.html
# https://vtk.org/doc/nightly/html/vtkDataWriter_8h_source.html
def convert_vtk_binary_to_ascii(input_file, output_file, fileformat="old"):
    print("vtk version = {}".format(vtk.vtkVersion.GetVTKVersion()))
    # Reader (auto-detects type: PolyData, UnstructuredGrid, etc.)
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(input_file)
    reader.Update()

    data = reader.GetOutput()

    # Writer
    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(data)
    writer.SetFileTypeToASCII()  # convert to ASCII
   
    print("write vtk version:{}".format(42))
    writer.SetFileVersion(42)  

    writer.Write()


if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    convert_vtk_binary_to_ascii(input_path, output_path)
