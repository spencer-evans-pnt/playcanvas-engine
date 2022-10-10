import { http } from '../platform/net/http.js';

class HtmlHandler {
    /**
     * Type of the resource the handler handles.
     *
     * @type {string}
     */
    handlerType = "html";

    constructor(app) {
        this.maxRetries = 0;
    }

    load(url, callback) {
        if (typeof url === 'string') {
            url = {
                load: url,
                original: url
            };
        }

        http.get(url.load, {
            retry: this.maxRetries > 0,
            maxRetries: this.maxRetries
        }, function (err, response) {
            if (!err) {
                callback(null, response);
            } else {
                callback(`Error loading html resource: ${url.original} [${err}]`);
            }
        });
    }

    open(url, data) {
        return data;
    }

    patch(asset, assets) {
    }
}

export { HtmlHandler };
