cwlVersion: v1.0
class: Workflow

inputs:
  input_file: File
  template_file: File

outputs:
  output_file:
    type: File
    outputSource: hello-world/output

steps:
  hello-world:
    run: dockstore-tool-helloworld.cwl
    in:
      input_file: input_file
      template_file: template_file
    out: [output]
