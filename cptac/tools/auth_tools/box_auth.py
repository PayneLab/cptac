#   Copyright 2018 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import logging
import requests
import threading
import webbrowser
from datetime import datetime, timedelta
from queue import Queue
from werkzeug import Request, Response
from werkzeug.serving import make_server

class BoxAuth:

    def __init__(self):
        # box token for accessing files via the Box API
        self.box_token = None
        
        # time delta for checking if a box token is still valid (they expire after an hour)
        self.box_token_age = None
        
        # a dictionary of passwords used for accessing password restricted files
        self.dataset_passwords = dict()
    
    def request_box_token(self):
        @Request.application
        def receive(request):
            q.put(request.args.get('code'))
            return Response("Authentication successful. You can close this window.", 200)

        # Don't show logs from server
        log = logging.getLogger('werkzeug')
        log.disabled = True

        # Set up authentication parameters
        base_url = "https://account.box.com/api/oauth2/authorize"
        client_id = "kztczhjoq3oes38yywuyfp4t9tu11it8"
        client_secret = "a5xNE1qj4Z4H3BSJEDVfzbxtmxID6iKY"
        login_url = f"{base_url}?client_id={client_id}&response_type=code"

        q = Queue()
        s = make_server("localhost", 8003, receive)
        t = threading.Thread(target=s.serve_forever)
        t.start()

        # Send the user to the "Grant access" page
        webbrowser.open(login_url)
        login_msg = "Please login to Box on the webpage that was just opened and grant access for cptac to download files through your account. If you accidentally closed the browser window, press Ctrl+C and call the download function again."
        print(login_msg)

        temp_code = q.get(block=True)
        s.shutdown()
        t.join()

        # Use the temporary access code to get the long term access token
        token_url = "https://api.box.com/oauth2/token";

        params = {
        'grant_type': 'authorization_code',
        'code': temp_code,
        'client_id': client_id,
        'client_secret': client_secret,
        }

        auth_resp = requests.post(token_url, data=params)
        self.box_token = auth_resp.json()["access_token"]
        self.box_token_age = datetime.now()

    def refresh_token(self):
        if self.box_token == None or datetime.now() - self.box_token_age > timedelta(minutes=59):
            self.request_box_token()

    def get_box_token(self):
        return self.box_token

    def get_box_token_age(self):
        return self.box_token_age

    def set_box_token(self, token, age):
        self.box_token = token
        self.box_token_age = age


    def add_single_password(self, dataset, password):
        self.dataset_passwords[dataset] = password
    
    def add_password_dictionary(self, passwords):
        if type(passwords) is dict:
            self.dataset_passwords = passwords

    def get_password(self, dataset):
        if dataset in self.dataset_passwords:
            return self.dataset_passwords[dataset]
        else:
            return None