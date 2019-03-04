# mmSCQC
Workflow for QC phase of Aaron Lun's simpleSingleCell Bioconductor workflow, expressed in WDL and CWL.

## Purpose

To demonstrate interface construction relating Bioconductor workflow packages/documents to workflows as construed by dockstore.org.

## Using cwltool

Clone the repository and verify that you have docker and cwltool installed.  The invocation
```
cwltool domm.cwl domm.yaml
```
should produce Rplots.pdf containing QC visualizations from the workflow.

## Using cromwell

```
java -jar cromwell-36.jar run mmSCQC.wdl
```
will produce, in the appropriate executions folder, plot1.pdf and plot2.pdf, the QC visualizations from the workflow.

## Other approaches

Through the dockstore.org API, the workflow can be executed locally using the dockstore CLI, or in hosted environments such as Firecloud, DNA Nexus, or DNAstack.  See the [dockstore workflow link](https://dockstore.org/workflows/github.com/vjcitn/mmSCQC/fromAaronLunSimpleSingleCell:master?tab=info); you will need to set up an account on the hosted platform to use this method.
