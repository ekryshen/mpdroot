# How to run gdb with ROOT macros

## Preparation
We are going to run gdb inside Visual Studio Code provided inside our toolbox.
One can (not mandatory) update Visual Studio code by command
```
sudo yum update code
```

## To start debugging (always)
1. enter c7-nica toolbox
1. `module add mpddev/v22.12.22-3`
1. `code`
1. Inside terminal of VS Code, source mpdroot file `source path_to_mpdroot/config/env.sh`

## Prepare json files (once)
settings.json
```
{
  "C_Cpp.default.includePath": [
    "${default}",
    "${env:ROOT_INCLUDE_PATH}"
  ],
  "files.associations": {
    "runMC.C": "cpp",
    "geometry_stage1.C": "cpp"
  },
}
```

launch.json
```
{
  "version": "0.2.0",
  "configurations": [
    {
      "name": "(gdb) Launch",
      "type": "cppdbg",
      "request": "launch",
      "program": "${env:ROOTSYS}/bin/root.exe",
      "args": [
        "-l",
        "-q",
        "${file}+g"
      ],
      "stopAtEntry": false,
      "cwd": "${workspaceFolder}",
      "environment": [],
      "externalConsole": false,
      "MIMode": "gdb",
      "setupCommands": [
        {
          "description": "Enable pretty-printing for gdb",
          "text": "-enable-pretty-printing",
          "ignoreFailures": true
        }
      ]
    }
  ]
}
```

## To debug
1. open macro file
1. place breakpoint
1. start debugging (via (gdb) Launch) - In the beginning launch will fail since header files in macro are missing. All missing header files need to be added as includes.

## Further reading
[View ROOT files inside VS Code](https://github.com/AlbertoPdRF/root-file-viewer#development)

[Run ROOT macros in VS Code](https://root.cern/blog/root-on-vscode/)
