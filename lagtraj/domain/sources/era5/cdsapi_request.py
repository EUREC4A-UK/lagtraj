import cdsapi
from cdsapi.api import Result


class RequestFetchCDSClient(cdsapi.Client):
    class RequestNotFoundException(Exception):
        pass

    """
    Wraps CDS api so that we can submit a request, get the request id and then
    later query the status or download data based on a request ID.
    """

    def __init__(self, *args, **kwargs):
        # LD: looking at the code forget=True avoids the cdsapi client sleeping
        # and waiting for the request to complete
        kwargs["forget"] = True
        super().__init__(*args, **kwargs)

    def queue_data_request(self, repository_name, query_kwargs):
        response = self.retrieve(repository_name, query_kwargs)

        if response.status_code not in [200, 202]:
            raise Exception(
                "Something went wrong requesting the data: {}"
                "".format(response.json())
            )
        else:
            reply = response.json()
            return reply["request_id"]

    def download_data_by_request(self, request_id, target):
        reply = self._get_request_status(request_id=request_id)

        result = Result(client=self, reply=reply)
        result.download(target=target)

    def _get_request_status(self, request_id):
        task_url = "{}/tasks/{}".format(self.url, request_id)
        session = self.session
        result = self.robust(session.get)(
            task_url, verify=self.verify, timeout=self.timeout
        )
        return result.json()

    def get_request_status(self, request_id):
        reply = self._get_request_status(request_id=request_id)
        if not "state" in reply:
            if reply["message"] == "Not found":
                raise self.RequestNotFoundException
            else:
                raise NotImplementedError(reply)
        else:
            return reply["state"]
