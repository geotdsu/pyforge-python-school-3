from redis import Redis
import json


redis_client = Redis(host="redis", port=6379, db=0)


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: list, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value))