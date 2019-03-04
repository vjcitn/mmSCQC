cwlVersion: v1.0
class: CommandLineTool
baseCommand: [R, --no-save]
hints:
  DockerRequirement:
    dockerPull: vjcitn/mmsc1
inputs:
  fcall:
    type: File?
    inputBinding:
      prefix: --file=
      separate: false
      position: 1
outputs:
  output:
    type:
      type: array
      items: File
    outputBinding: 
      glob: '*.pdf'
