# How to Download and Run the Executable

## 1. Download the Release

- Go to the [Releases](../../releases) section of this repository.
- Download the latest release package (e.g., `YourApp_v1.0.zip`).
- Extract the zip file to a folder of your choice.

## 2. Install the MATLAB Runtime

This executable requires the free [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html) (MCR) to be installed.

**Steps:**
1. Visit the [MATLAB Runtime download page](https://www.mathworks.com/products/compiler/matlab-runtime.html).
2. Download the version **matching the MATLAB version used to compile this app** (e.g., R2024a).  
   *The required version is usually mentioned in the release notes or by the developer.*
3. Run the installer and follow the on-screen instructions to complete the installation.

## 3. Run the Executable

- Navigate to the extracted release folder.
- Double-click the executable file (e.g., `YourApp.exe` on Windows, `YourApp` on Linux/Mac).
- If you see a message about missing MATLAB Runtime, ensure you have installed the correct version.

### Command Line (optional)

You can also run the executable from the command line:
```sh
# Windows
YourApp.exe

# Linux/Mac
./YourApp
```

## 4. Troubleshooting

- **Error: Missing MATLAB Runtime**  
  Make sure you have installed the correct version of MATLAB Runtime.
- **Application does not start**  
  Check that your system meets the minimum requirements and that you have extracted all files from the zip archive.

## 5. References

- [MATLAB Runtime FAQ](https://www.mathworks.com/help/compiler/matlab-runtime.html)
- [MATLAB Compiler Documentation](https://www.mathworks.com/help/compiler/)
