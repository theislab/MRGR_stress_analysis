# Single-cell characterization of the stress response in the adult ventral hippocampus

> scRNA-seq analyses

[![](/readme.png)]()

## Structure
```
├── README.md
├── docker
│   ├── Dockerfile
│   └── python-packages.txt
├── notebooks
│   ├── 01-QC
│   │   ├── 01a-QC_GrCtrl.ipynb
│   │   ├── 01b-QC_GrDlx.ipynb
│   │   ├── 01c-QC_GrNex.ipynb
│   │   ├── 01d-QC_MrCtrl.ipynb
│   │   ├── 01e-QC_MrDlx.ipynb
│   │   ├── 01f-QC_MrNex.ipynb
│   │   └── 02-QC_DoubletScran.ipynb
│   ├── 02-Annotation
│   │   ├── 01-Annotation.ipynb
│   │   └── 02-AnatomicalAnnotation.ipynb
│   ├── 03-DE
│   │   ├── 01a-DEPostFINAL.ipynb
│   │   └── helper_de.py
│   ├── 05-Trajectory
│   │   └── 01-CellrankDPT.ipynb
│   └── 06-Paper
│       ├── Figure1.ipynb
│       ├── Figure2-3-4.ipynb
│       ├── Figure5-upset.ipynb
│       ├── Figure5.ipynb
│       ├── Figure6.ipynb
│       └── Figure7.ipynb
└── scripts
    └── 00_diffxpy_stress.py
```

## Installation

We provide a Dockerfile with all the packages required to run the analyses. To build the Docker image use the command:
```shell
docker build docker -t mrgr:final
```

After the image is compiled, you can run it interactively using:

```shell
docker run --interactive --tty --name mrgr --publish 8888:8888 --volume $HOME:/root/host_home --workdir /root mrgr:final  /bin/bash
```

This will start a container, in order to start a JupyterLab session use the alias:

```shell
jl
```

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**