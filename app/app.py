import os
import json
import pathlib
this_dir = pathlib.Path(__file__).parent.resolve()

from flask import Flask
app = Flask(__name__)



@app.route('/hello')
def hello_world():
   return 'Hello World flask'

if __name__ == '__main__':
   

   with open(os.path.join(this_dir, '.config_cache.json'), 'r') as f:
      config = json.load(f)

   app.run(port=int(config['flask_port']),debug=True)
