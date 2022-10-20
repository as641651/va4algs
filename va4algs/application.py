import json
import os

class Application:
    def __init__(self, 
        local_problem_dir,
        flask_port,
        default_operand_size,
        language,
        app_dir="app/"):

        self.local_problem_dir = local_problem_dir
        self.flask_port = flask_port
        self.default_operand_size = default_operand_size
        self.language = language
        self.app_dir = app_dir
        self.m_server = -1
        self.m_uname = -1
        self.m_problem_dir = -1
        self.m_init_script = -1

    def use_measurement_server(self, server, 
        uname, 
        problem_dir,
        init_script):
        
        self.m_server = server
        self.m_uname = uname
        self.m_problem_dir = problem_dir
        self.m_init_script = init_script 

    def initialize(self):
        config = {
            "local_problem_dir":self.local_problem_dir,
            "flask_port":self.flask_port,
            "default_operand_size": self.default_operand_size,
            "language":self.language,
            "m_server":self.m_server,
            "m_uname":self.m_uname,
            "m_problem_dir":self.m_problem_dir,
            "m_init_script":self.m_init_script
        }

        json_object = json.dumps(config, indent=4)
        with open(os.path.join(self.app_dir,".config_cache.json"), "w") as outfile:
            outfile.write(json_object)


