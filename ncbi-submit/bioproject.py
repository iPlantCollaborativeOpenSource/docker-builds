__author__ = 'Dennis Roberts'

import config.ncbi_submit_properties
from lxml import etree

def load_schema(schemas_dir, path):
    with open(path, 'r') as f:
        schema_root = etree.XML(f.read(), base_url=schemas_dir)
    return etree.XMLSchema(schema_root);

def load_parser(schemas_dir, path):
    schema = load_schema(schemas_dir, path)
    return etree.XMLParser(schema=schema)

class BioProjectXmlValidator:
    def __init__(self, schemas_dir, schema_paths):
        self.xmlparser = load_parser(schemas_dir, schema_paths['submission'])
        self.bioproject_schema = load_schema(schemas_dir, schema_paths['bioproject'])
        self.biosample_schema = load_schema(schemas_dir, schema_paths['biosample'])
        self.genome_schema = load_schema(schemas_dir, schema_paths['genome'])

    def validate_bioproject_xml(self, submission_path):
        with open(submission_path, 'r') as f:
            submission = etree.fromstring(f.read(), self.xmlparser)

        bioproject = submission.xpath('/Submission/Action/AddData/Data/XmlContent/Project')
        if len(bioproject) > 0:
            self.bioproject_schema.assertValid(bioproject[0])

        bio_samples = submission.xpath('/Submission/Action/AddData/Data/XmlContent/BioSample')
        for biosample in bio_samples:
            self.biosample_schema.assertValid(biosample)

class BioProjectUploader:
    def __init__(self, ascp_cmd, private_key_path, ncbi_user, ncbi_host, ncbi_sumbit_path):
        self.ascp_cmd = ascp_cmd
        self.private_key_path = private_key_path
        self.upload_dest = '{0}@{1}:{2}'.format(ncbi_user, ncbi_host, ncbi_sumbit_path)

    def upload_project(self, submit_dir, input_paths):
        # Collect the submission files from the input paths into the submission dir
        src_files = {}
        for path in input_paths:
            filename = os.path.basename(path)
            if filename in src_files:
                fmt = "Duplicate filenames found in input directory:\n{0}\n{1}"
                raise Exception(fmt.format(src_files[filename], path))
            src_files[filename] = path

            shutil.move(path, os.path.join(submit_dir, filename))

        ascp_cmd = self.ascp_cmd + [
            "-i", self.private_key_path,
            submit_dir,
            self.upload_dest
        ]

        try:
            retcode = call(ascp_cmd)
            if retcode != 0:
                raise Exception("Upload error: {0}".format(-retcode))

            # The file uploads were successful, so upload a 'submit.ready' file to complete the submission.
            submit_ready = "submit.ready"
            open(os.path.join(submit_dir, submit_ready), 'a').close()

            # Calling the same upload command with the same submit directory will skip all files already
            # successfully uploaded, and only upload the new 'submit.ready' file.
            retcode = call(ascp_cmd)
            if retcode != 0:
                raise Exception("Error uploading '{0}' file: {1}".format(submit_ready, -retcode))
        except OSError as e:
            raise Exception("Aspera Connect upload failed", e)

        # Clean up: Move input files back into their original directories,
        # so they are not transferred as outputs, but can be preserved as inputs
        for filename in src_files:
            shutil.move(os.path.join(submit_dir, filename), src_files[filename])

def get_xml_validator():
    return BioProjectXmlValidator(
        config.ncbi_submit_properties.schemas_dir,
        config.ncbi_submit_properties.schema_paths
    )

def get_uploader(private_key_path=None):
    key_path = config.ncbi_submit_properties.private_key_path if private_key_path is None else private_key_path
    BioProjectUploader(
        config.ncbi_submit_properties.ascp_cmd,
        key_path,
        config.ncbi_submit_properties.ncbi_user,
        config.ncbi_submit_properties.ncbi_host,
        config.ncbi_submit_properties.ncbi_sumbit_path
    )