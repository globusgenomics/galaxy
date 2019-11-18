# Globus Genomics Galaxy

S3 + Kubernetes setup:
Galaxy's job execution requires job working directory, inputs, outputs, tools and indices available on the local filesystem. In order to use distributed computation resources, NFS has to be mounted to the workers so that these files are accessible. But NFS has its limitation, it creates large network burden for the head node when there are a lot of workers running at the same time. We've developed a solution to solve this issue and in the meantime lower the cost.

S3 is a low cost storage, it is durable and highly available, also data transfer is free and fast inside the same AWS region, we only have to pay the data operation cost which is reasonable. So S3 is very ideal to be used as the main storage and a relay to distribute data among the workers. Also a S3 bucket can be mounted as file system, but this is only stable enough to be used as read only file system, so we won't use it for workers, only for the head node when users want to read output datasets.

Containers keep a clean environment for the job executions. Kubernetes is a great tool to help with running containers on the cloud. It comes with a great auto scaler to scale up and down the worker fleets. The auto scaler can launch workers based on job requirements and also diversify the instance types. This increases the resource utilisation and lower the chance of job evictions, since we use spot instances to lower the cost and having diversified instance types is very helpful.

Development: Galaxy's Kubernetes runner is updated and a job execution script is created to achieve this goal. 
- Pre job: On the head node, from the job execution files, get the job working directory path and sync the directory to S3; Get the inputs and outputs paths, upload them to S3 if not already exist on S3, create symlinks pointing to the file in the mounted S3 file system; Get the tools information and indices information; Replace the job execution command with the command to run the job execution script along with the arguments.
- Run job: In the container on the worker, the job execution script will download all the required working directory, inputs and outputs from S3; Tools are provided by NFS; Indices will also be downloaded from S3; Execute the job; After the job is finished, upload outputs and working directory back to S3, clean out the data.
- Post job: On the head node, sync the job working directory from S3 so that the job logs and exit code files are available on the local file system. The users can read the outputs because the S3 is mounted to the head node.

To use this setup, the tool needs to be configured to use the Kubernetes runner in the job_conf.xml, that is choose one of the k8s destinations depends on the job resource requirements; Provide S3 bucket information in the galaxy.yml, and make sure head node and workers have access to the bucket such as using IAM roles.

