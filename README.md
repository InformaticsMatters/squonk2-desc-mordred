# squonk2-desc-mordred

A Squonk2 Data Manager Job that calculates [Mordred] descriptors for a set of
molecules, using the [mordredcommunity] fork and RDKit.

The Job Definition lives in [`data-manager/jobs.yaml`](data-manager/jobs.yaml);
documentation for the Job is in [`data-manager/docs`](data-manager/docs).

## The image

The image is built from [`Dockerfile`](Dockerfile) and published as
`informaticsmatters/mordred-descriptors`. Dependencies are managed with Poetry;
`poetry.lock` is what the image installs, so a dependency change means
relocking with the Poetry version the `Dockerfile` pins.

## Testing

Jobs are tested with [jote]:

```bash
poetry install --only dev
jote
```

`jote --dry-run` validates the Job Definitions against the schema without
running anything, which is what CI does on every branch.

## Releasing

1. Run the `publish-tag` workflow with the new tag.
2. Set the Job Definition `version` and `image.tag` to that same value.
3. Cut the matching Git tag.

Never reuse a container tag — the Data Manager caches any tag other than
`latest`/`stable` per Kubernetes node. See `docs/versioning.md` in the
[squonk2-jobs] umbrella repository.

[jote]: https://github.com/InformaticsMatters/squonk2-data-manager-job-tester
[mordred]: https://github.com/mordred-descriptor/mordred
[mordredcommunity]: https://github.com/JacksonBurns/mordred-community
[squonk2-jobs]: https://github.com/InformaticsMatters/squonk2-jobs
